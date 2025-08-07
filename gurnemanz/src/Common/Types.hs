{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE OverloadedStrings #-}

module Common.Types
    ( -- * Resource Types
      LogHandle(..)
    , Resources(..)
    -- * Type Classes
    , HasResources(..)
    , HasConfig(..)
    -- * Log Types
    , LogLevel(..)
    , LogFormat(..)
    , LogConfig(..)
    , LogContext(..)
    -- * Config Types
    , AppConfig(..)
    , Environment(..)
    , ChemConfig(..)
    , APIConfig(..)
    -- * Type Re-exports
    , Word16
    ) where

import Control.Concurrent (MVar)
import Data.Text (Text)
import Data.Word (Word16)
import Data.Aeson (FromJSON(..), ToJSON(..), withText, Value(String))
import GHC.Generics (Generic)
import System.IO (Handle)


-- | Log levels for the application
data LogLevel
    = Debug
    | Info
    | Warning
    | Error
    deriving stock (Show, Eq, Ord, Generic)
    deriving anyclass (ToJSON, FromJSON)

-- | Log format specification
data LogFormat
    = LogText      -- ^ Plain text logs
    | LogJSON      -- ^ JSON formatted logs
    deriving stock (Show, Eq, Generic)

instance ToJSON LogFormat where
    toJSON LogText = String "text"
    toJSON LogJSON = String "json"

instance FromJSON LogFormat where
    parseJSON = withText "LogFormat" $ \case
        "text" -> pure LogText
        "json" -> pure LogJSON
        _ -> fail "Invalid log format"

-- | Logging context
data LogContext = LogContext
    { contextLevel :: !LogLevel
    , contextComponent :: !Text
    , contextCorrelationId :: !(Maybe Text)
    , contextTags :: ![(Text, Text)]
    } deriving stock (Show, Eq, Generic)
      deriving anyclass (ToJSON, FromJSON)

-- | Logging configuration
data LogConfig = LogConfig
    { logLevel :: !LogLevel
    , logFile :: !(Maybe FilePath)
    , logFormat :: !LogFormat
    , enableConsoleLog :: !Bool
    } deriving stock (Show, Eq, Generic)
      deriving anyclass (FromJSON, ToJSON)

instance Semigroup LogConfig where
    a <> b = LogConfig
        { logLevel = max (logLevel a) (logLevel b)
        , logFile = logFile b <|> logFile a
        , logFormat = logFormat b
        , enableConsoleLog = enableConsoleLog a || enableConsoleLog b
        }

instance Monoid LogConfig where
    mempty = LogConfig
        { logLevel = Info
        , logFile = Nothing
        , logFormat = LogText
        , enableConsoleLog = False
        }

(<|>) :: Maybe a -> Maybe a -> Maybe a
Nothing <|> b = b
a <|> _ = a

-- | Logging handle wrapper
newtype LogHandle = LogHandle
    { unLogHandle :: MVar Handle }


-- | Chemistry configuration
data ChemConfig = ChemConfig
    { maxMoleculeSize :: !Int        -- ^ Maximum size of molecules
    , maxTransformationSteps :: !Int  -- ^ Maximum transformation chain length
    , computationTimeout :: !Int      -- ^ Computation timeout in seconds
    , enabledProperties :: ![Text]    -- ^ List of enabled property calculations
    } deriving stock (Show, Eq, Generic)
      deriving anyclass (FromJSON, ToJSON)

-- | API configuration
data APIConfig = APIConfig
    { apiPort :: !Word16       -- ^ API server port
    , apiHost :: !Text         -- ^ API server host
    , apiMaxRequestSize :: !Int  -- ^ Maximum request size in bytes
    , apiRateLimit :: !Int      -- ^ Rate limit per minute
    , apiTimeout :: !Int        -- ^ Request timeout in seconds
    , socketPath :: !FilePath   -- ^ Unix domain socket path
    } deriving stock (Show, Eq, Generic)
      deriving anyclass (FromJSON, ToJSON)

-- | Application environment
data Environment
    = Development   -- ^ Development environment
    | Staging      -- ^ Staging environment
    | Production   -- ^ Production environment
    deriving stock (Show, Eq, Generic)
    deriving anyclass (FromJSON, ToJSON)

-- | Application configuration
data AppConfig = AppConfig
    { chemConfig :: !ChemConfig      -- ^ Chemical computation configuration
    , apiConfig :: !APIConfig        -- ^ API server configuration
    , logConfig :: !LogConfig        -- ^ Logging configuration
    , environment :: !Environment    -- ^ Application environment
    } deriving stock (Show, Eq, Generic)
      deriving anyclass (FromJSON, ToJSON)

-- | Application resources
data Resources = Resources
    { resourceLogHandle :: !LogHandle               -- ^ Logging handle
    }

-- | Access to resources
class HasResources env where
    getResources :: env -> Resources

-- | Access to configuration
class HasConfig env where
    getConfig :: env -> AppConfig