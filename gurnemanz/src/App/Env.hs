{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}

module App.Env
    ( -- * Environment Management
      initializeEnv
    , shutdownEnv
    , withAppEnv
    ) where

import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Error.Class (MonadError, throwError)
import Control.Exception (try, SomeException, bracket, throwIO, toException)
import Control.Monad.Except (runExceptT)
import Data.Text (Text)
import qualified Data.Text as T


import Core.Error (Error(..))
import Infrastructure.Config.Types (AppConfig(..))
import Common.Types (LogConfig(..), LogLevel(..), LogContext(..), Environment(..))
import Infrastructure.Config.Environment (validateEnvironment)
import Core.Types (Resources(..))
import App.Types (AppEnv(..))
import qualified Infrastructure.Persistence.Logging as Log

-- | Initialize application environment
initializeEnv :: (MonadIO m, MonadError Error m) => AppConfig -> m AppEnv
initializeEnv config = do
    -- Validate environment configuration
    validateEnvironment config

    -- Initialize resources
    resources <- initializeResources config

    pure AppEnv
        { appConfig = config
        , appResources = resources
        }

-- | Shutdown application environment
shutdownEnv :: MonadIO m => AppEnv -> m ()
shutdownEnv = shutdownResources

-- | Bracket pattern for environment management
withAppEnv :: AppConfig -> (AppEnv -> IO a) -> IO a
withAppEnv config action = bracket
    (runExceptT $ initializeEnv config)
    (either (const $ pure ()) (liftIO . shutdownEnv))
    (either (throwIO . toException) action)

-- | Initialize all resources
initializeResources :: (MonadIO m, MonadError Error m) => AppConfig -> m Resources
initializeResources AppConfig{..} = do
    -- Initialize logging
    logHandle <- Log.initializeLogging logConfig

    pure Resources
        { resourceLogHandle = logHandle
        }


-- | Shutdown all resources
shutdownResources :: MonadIO m => AppEnv -> m ()
shutdownResources AppEnv{..} = do
    -- Shutdown logging
    Log.shutdownLogging (resourceLogHandle appResources)