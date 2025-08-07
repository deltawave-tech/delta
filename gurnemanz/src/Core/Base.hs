{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleInstances #-}

module Core.Base
    ( -- * Type Aliases
      UserId
    , SessionId
    , MolId
    , MoleculeId
    , TransformationId
    , MethodId(..)
    , Result
    -- * Logging Types (re-exported)
    , LogLevel(..)
    , LogFormat(..)
    , LogConfig(..)
    , LogContext(..)
    -- * Logger Interface
    , MonadLogger(..)
    ) where

import Data.Text (Text)
import Data.UUID (UUID)
import Control.Monad.IO.Class (MonadIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Reader (ReaderT)
import Control.Monad.Except (ExceptT)
import Data.Aeson (ToJSON(..), FromJSON(..), withText, Value(String), object)
import GHC.Generics (Generic)

-- Import from Common.Types to avoid circular dependencies
import Common.Types 
    ( LogLevel(..)
    , LogFormat(..)
    , LogConfig(..)
    , LogContext(..)
    )

-- | Basic type aliases
type UserId = Text
type SessionId = UUID
type MolId = UUID
type MoleculeId = UUID
type TransformationId = UUID
type Result = ()

-- | Method identifier type
data MethodId
    = Method Text
    deriving stock (Show, Eq, Generic)
    deriving anyclass (ToJSON, FromJSON)

-- | Logger interface
class MonadIO m => MonadLogger m where
    logInfo :: Text -> m ()
    logWarning :: Text -> m ()
    logError :: Text -> m ()
    logDebug :: Text -> m ()

instance MonadLogger m => MonadLogger (ExceptT e m) where
    logInfo = lift . logInfo
    logWarning = lift . logWarning
    logError = lift . logError
    logDebug = lift . logDebug

instance MonadLogger m => MonadLogger (ReaderT r m) where
    logInfo = lift . logInfo
    logWarning = lift . logWarning
    logError = lift . logError
    logDebug = lift . logDebug