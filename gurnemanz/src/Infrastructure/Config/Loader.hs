{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}

module Infrastructure.Config.Loader
    ( -- * Configuration Loading
      loadConfig
    , validateConfig
      -- * Default Configurations
    , defaultAppConfig
    , defaultChemConfig
    , defaultAPIConfig
    , defaultLogConfig
      -- * Configuration Validation
    , validateChemConfig
    , validateLogConfig
    , validateAPIConfig
    ) where

import qualified Data.Yaml as Y
import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Error.Class (MonadError)
import Data.Text (Text)
import qualified Data.Text as T
import Control.Monad (when, unless)
import Data.Word (Word16)

import Core.Error (Error(..), throwError)
import Infrastructure.Config.Types
    ( AppConfig(..)
    , Environment(..)
    , ChemConfig(..)
    , APIConfig(..)
    )
import Core.Base (LogLevel(..), LogFormat(..), LogConfig(..))


-- | Default chemical computation configuration
defaultChemConfig :: ChemConfig
defaultChemConfig = ChemConfig
    { maxMoleculeSize = 1000
    , maxTransformationSteps = 50
    , computationTimeout = 300
    , enabledProperties =
        [ "logP"
        , "topology"
        , "molecular_weight"
        , "rotatable_bonds"
        ]
    }

-- | Default API configuration
defaultAPIConfig :: APIConfig
defaultAPIConfig = APIConfig
    { apiPort = 8080
    , apiHost = "0.0.0.0"
    , apiMaxRequestSize = 10 * 1024 * 1024  -- 10MB
    , apiRateLimit = 100  -- requests per minute
    , apiTimeout = 30     -- seconds
    , socketPath = "/tmp/gurnemanz.sock"
    }

-- | Default logging configuration
defaultLogConfig :: LogConfig
defaultLogConfig = LogConfig
    { logLevel = Info
    , logFile = Just "app.log"
    , logFormat = LogText
    , enableConsoleLog = True
    }

-- | Default complete application configuration
defaultAppConfig :: AppConfig
defaultAppConfig = AppConfig
    { chemConfig = defaultChemConfig
    , apiConfig = defaultAPIConfig
    , logConfig = defaultLogConfig
    , environment = Development
    }

-- | Load configuration from file
loadConfig :: (MonadIO m, MonadError Error m) => FilePath -> m AppConfig
loadConfig path = do
    contents <- liftIO $ Y.decodeFileEither path
    case contents of
        Left err -> throwError $ ValidationError $
            "Could not parse config: " <> T.pack (show err)
        Right config -> validateConfig config

-- | Validate complete configuration
validateConfig :: MonadError Error m => AppConfig -> m AppConfig
validateConfig config@AppConfig{..} = do
    validateChemConfig chemConfig
    validateLogConfig logConfig
    validateAPIConfig apiConfig
    pure config


-- | Validate chemical computation configuration
validateChemConfig :: MonadError Error m => ChemConfig -> m ()
validateChemConfig ChemConfig{..} = do
    when (maxMoleculeSize <= 0) $
        throwError $ ValidationError "Maximum molecule size must be positive"
    when (maxTransformationSteps <= 0) $
        throwError $ ValidationError "Maximum transformation steps must be positive"
    when (computationTimeout <= 0) $
        throwError $ ValidationError "Computation timeout must be positive"
    unless (not $ null enabledProperties) $
        throwError $ ValidationError "At least one property must be enabled"

-- | Validate logging configuration
validateLogConfig :: MonadError Error m => LogConfig -> m ()
validateLogConfig LogConfig{..} = do
    case logFile of
        Just path -> when (null path) $
            throwError $ ValidationError "Log file path cannot be empty"
        Nothing -> when (not enableConsoleLog) $
            throwError $ ValidationError "Must enable either file or console logging"

-- | Validate API configuration
validateAPIConfig :: MonadError Error m => APIConfig -> m ()
validateAPIConfig APIConfig{..} = do
    validatePort apiPort "API"
    when (apiMaxRequestSize <= 0) $
        throwError $ ValidationError "API max request size must be positive"
    when (apiRateLimit <= 0) $
        throwError $ ValidationError "API rate limit must be positive"
    when (apiTimeout <= 0) $
        throwError $ ValidationError "API timeout must be positive"
    when (null socketPath) $
        throwError $ ValidationError "Socket path cannot be empty"

-- | Validate port number
validatePort :: MonadError Error m => Word16 -> Text -> m ()
validatePort port service =
    when (port <= 0 || port > 65535) $
        throwError $ ValidationError $
            service <> " port must be between 1 and 65535"
