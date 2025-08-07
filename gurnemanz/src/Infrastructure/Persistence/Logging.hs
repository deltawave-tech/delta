{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE UndecidableInstances #-}

module Infrastructure.Persistence.Logging
    ( -- * Logging initialization
      initializeLogging
    , shutdownLogging
      -- * Logging operations
    , logMessage
      -- * Logging types
    , LogTarget(..)
    , LogHandle(..)
      -- * Type classes
    , HasLogHandle(..)
    , HasLogConfig(..)
    , HasLogContext(..)
      -- * Logging monad
    , LoggingT(..)
    , runLogging
    ) where

import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Reader (MonadReader, local)
import Control.Monad.Error.Class (MonadError(..))
import Control.Concurrent (MVar, newMVar, readMVar, withMVar)
import Data.Text (Text)
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import Data.Time (UTCTime, getCurrentTime, formatTime, defaultTimeLocale)
import System.IO (Handle, IOMode(..), openFile, hClose, stdout)
import qualified System.IO as IO
import Data.Aeson (ToJSON(..), FromJSON(..), encode)
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Text.Encoding (decodeUtf8)
import Control.Monad (unless, when)
import GHC.Generics (Generic)

import Core.Error (Error(..))
import Common.Types
    ( LogLevel(..)
    , LogConfig(..)
    , LogFormat(..)
    , LogContext(..)
    , LogHandle(..)
    )

-- | Basic logging monad
newtype LoggingT m a = LoggingT
    { runLoggingT :: LogHandle -> LogConfig -> LogContext -> m a }

-- Manual instances
instance Functor m => Functor (LoggingT m) where
    fmap f (LoggingT r) = LoggingT $ \h c ctx ->
        fmap f (r h c ctx)

instance Applicative m => Applicative (LoggingT m) where
    pure a = LoggingT $ \_ _ _ -> pure a
    LoggingT f <*> LoggingT a = LoggingT $ \h c ctx ->
        f h c ctx <*> a h c ctx

instance Monad m => Monad (LoggingT m) where
    LoggingT ma >>= f = LoggingT $ \h c ctx -> do
        a <- ma h c ctx
        runLoggingT (f a) h c ctx

instance MonadIO m => MonadIO (LoggingT m) where
    liftIO ma = LoggingT $ \_ _ _ -> liftIO ma

instance MonadError e m => MonadError e (LoggingT m) where
    throwError e = LoggingT $ \_ _ _ -> throwError e
    catchError (LoggingT m) h = LoggingT $ \handle config ctx ->
        catchError (m handle config ctx) (\e ->
            runLoggingT (h e) handle config ctx)

-- | Run logging operations
runLogging :: LogHandle -> LogConfig -> LogContext -> LoggingT IO a -> IO a
runLogging handle config ctx action = runLoggingT action handle config ctx

-- | Logging target configuration
data LogTarget
    = LogFile FilePath
    | LogStdout
    | LogBoth FilePath
    deriving stock (Show, Eq, Generic)
    deriving anyclass (ToJSON, FromJSON)

-- | Log entry structure
data LogEntry = LogEntry
    { entryTimestamp :: !UTCTime
    , entryLevel :: !LogLevel
    , entryMessage :: !Text
    , entryContext :: !LogContext
    , entryError :: !(Maybe Error)
    } deriving stock (Show, Generic)
      deriving anyclass (ToJSON)

-- | Initialize logging system
initializeLogging :: (MonadIO m, MonadError Error m) => LogConfig -> m LogHandle
initializeLogging LogConfig{..} = liftIO $ do
    handle <- case logFile of
        Nothing -> pure stdout
        Just path -> openFile path AppendMode
    unless (handle == stdout) $
        IO.hSetBuffering handle IO.LineBuffering
    LogHandle <$> newMVar handle

-- | Shutdown logging system
shutdownLogging :: MonadIO m => LogHandle -> m ()
shutdownLogging (LogHandle handleVar) = liftIO $ do
    handle <- readMVar handleVar
    unless (handle == stdout) $
        hClose handle

-- | Log a message with context
logMessage
    :: MonadIO m
    => LogLevel
    -> Text
    -> Maybe Error
    -> LogHandle
    -> LogConfig
    -> LogContext
    -> m ()
logMessage level msg err handle config ctx =
    when (level >= logLevel config) $ liftIO $ do
        timestamp <- getCurrentTime
        let entry = LogEntry
                { entryTimestamp = timestamp
                , entryLevel = level
                , entryMessage = msg
                , entryContext = ctx
                , entryError = err
                }
        withMVar (unLogHandle handle) $ \h ->
            writeLogEntry h (logFormat config) entry

-- | Format log entry as text
formatLogText :: LogEntry -> Text
formatLogText LogEntry{..} = T.concat
    [ formatTimestamp entryTimestamp
    , " [" <> formatLevel entryLevel <> "] "
    , entryMessage
    , maybe "" ((" | " <>) . T.pack . show) entryError
    ]

-- | Format log entry as JSON
formatLogJSON :: LogEntry -> Text
formatLogJSON entry =
    decodeUtf8 $ BL.toStrict $ encode entry

-- Helper functions
writeLogEntry :: Handle -> LogFormat -> LogEntry -> IO ()
writeLogEntry handle format entry = do
    let text = case format of
            LogText -> formatLogText entry
            LogJSON -> formatLogJSON entry
    TIO.hPutStrLn handle text

formatTimestamp :: UTCTime -> Text
formatTimestamp = T.pack . formatTime defaultTimeLocale "%Y-%m-%d %H:%M:%S"

formatLevel :: LogLevel -> Text
formatLevel = T.pack . show

-- Type classes for configuration access
class HasLogHandle a where
    getLogHandle :: a -> LogHandle

class HasLogConfig a where
    getLogConfig :: a -> LogConfig

class HasLogContext a where
    getLogContext :: a -> LogContext
    setLogContext :: LogContext -> a -> a

-- Context management
withLogContext
    :: (MonadReader r m, HasLogContext r)
    => LogContext
    -> m a
    -> m a
withLogContext newCtx = local (setLogContext newCtx)

modifyLogContext
    :: (MonadReader r m, HasLogContext r)
    => (LogContext -> LogContext)
    -> m a
    -> m a
modifyLogContext f = local (\r -> setLogContext (f (getLogContext r)) r)