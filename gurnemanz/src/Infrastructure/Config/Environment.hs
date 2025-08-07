{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}

module Infrastructure.Config.Environment
    ( -- * Environment Helpers
      isDevelopment
    , isProduction
    , isStaging
      -- * Environment-specific configs
    , getEnvironmentConfig
    , defaultDevelopmentConfig
    , defaultProductionConfig
    , defaultStagingConfig
      -- * Environment validation
    , validateEnvironment
    ) where

import Control.Monad.Error.Class (MonadError)
import Control.Monad (when)

import Core.Error (Error(..), throwError)
import Core.Base (LogLevel(..), LogConfig(..))
import Infrastructure.Config.Types
    ( AppConfig(..)
    , Environment(..)
    )
import Infrastructure.Config.Loader (defaultAppConfig, defaultLogConfig)

-- | Check if environment is Development
isDevelopment :: AppConfig -> Bool
isDevelopment config = environment config == Development

-- | Check if environment is Production
isProduction :: AppConfig -> Bool
isProduction config = environment config == Production

-- | Check if environment is Staging
isStaging :: AppConfig -> Bool
isStaging config = environment config == Staging

-- | Get environment-specific configuration
getEnvironmentConfig :: Environment -> AppConfig
getEnvironmentConfig = \case
    Development -> defaultDevelopmentConfig
    Production -> defaultProductionConfig
    Staging -> defaultStagingConfig

-- | Default development configuration
defaultDevelopmentConfig :: AppConfig
defaultDevelopmentConfig = defaultAppConfig
    { environment = Development
    , logConfig = (defaultLogConfig)
        { logLevel = Debug
        , enableConsoleLog = True
        }
    }

-- | Default production configuration
defaultProductionConfig :: AppConfig
defaultProductionConfig = defaultAppConfig
    { environment = Production
    , logConfig = (defaultLogConfig)
        { logLevel = Warning
        , enableConsoleLog = False
        }
    }

-- | Default staging configuration
defaultStagingConfig :: AppConfig
defaultStagingConfig = defaultAppConfig
    { environment = Staging
    , logConfig = (defaultLogConfig)
        { logLevel = Info
        , enableConsoleLog = True
        }
    }

-- | Validate environment-specific configuration
validateEnvironment :: MonadError Error m => AppConfig -> m ()
validateEnvironment AppConfig{..} = do
    -- Validate production settings
    when (isProduction AppConfig{..}) $ do
        when (logLevel logConfig < Warning) $
            throwError $ ValidationError "Production environment requires at least Warning log level"
        when (enableConsoleLog logConfig) $
            throwError $ ValidationError "Production environment should not enable console logging"

    -- Validate development settings
    when (isDevelopment AppConfig{..}) $ do
        when (not $ enableConsoleLog logConfig) $
            throwError $ ValidationError "Development environment should enable console logging"

    -- Validate staging settings
    when (isStaging AppConfig{..}) $ do
        when (logLevel logConfig > Debug) $
            throwError $ ValidationError "Staging environment should enable detailed logging"
