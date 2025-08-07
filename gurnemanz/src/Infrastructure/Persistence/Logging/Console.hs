{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE CPP #-}

module Infrastructure.Persistence.Logging.Console 
    ( withConsoleLogger
    , logToConsole 
    ) where

import Control.Monad.IO.Class (MonadIO, liftIO)
import Data.Time (getCurrentTime, formatTime, defaultTimeLocale)
import System.IO (stdout, hFlush)
import qualified Data.Text as T
import qualified Data.Text.IO as TIO

import Core.Base (MonadLogger(..))

-- | Console logging for IO computations
instance MonadLogger IO where
    logDebug msg = logToConsole "DEBUG" msg
    logInfo msg = logToConsole "INFO" msg
    logWarning msg = logToConsole "WARNING" msg
    logError msg = logToConsole "ERROR" msg

-- | Log directly to the console
logToConsole :: MonadIO m => String -> T.Text -> m ()
logToConsole level msg = liftIO $ do
    now <- getCurrentTime
    let timestamp = T.pack $ formatTime defaultTimeLocale "%Y-%m-%d %H:%M:%S" now
    let prefix = T.pack $ "[" ++ level ++ "]"
    TIO.putStrLn $ timestamp <> " " <> prefix <> " " <> msg
    hFlush stdout
    
-- | Run an action with console logging enabled
withConsoleLogger :: IO a -> IO a
withConsoleLogger action = do
    logToConsole "INFO" "Console logger initialized"
    result <- action
    logToConsole "INFO" "Console logger shutting down"
    return result