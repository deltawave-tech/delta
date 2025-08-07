{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveAnyClass #-}

module Infrastructure.Logging.Types
    ( LogEntry(..)
    , LogFormat(..)
    , LogLevel(..)
    , LogContext(..)
    ) where

import Data.Text (Text)
import Data.Time (UTCTime)
import Data.Aeson (ToJSON(..), FromJSON(..))
import GHC.Generics (Generic)
import Core.Error (Error)

-- | Log levels
data LogLevel
    = Debug     -- ^ Debug level logging
    | Info      -- ^ Information level logging
    | Warning   -- ^ Warning level logging
    | Error     -- ^ Error level logging
    deriving stock (Show, Eq, Ord, Generic)
    deriving anyclass (FromJSON, ToJSON)

-- | Logging context
data LogContext = LogContext
    { contextLevel :: !LogLevel
    , contextComponent :: !Text
    , contextCorrelationId :: !(Maybe Text)
    , contextTags :: ![(Text, Text)]
    } deriving stock (Show, Eq, Generic)
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

-- | Log format specification
data LogFormat
    = LogText      -- ^ Plain text logs
    | LogJSON      -- ^ JSON formatted logs
    | LogCustom (LogEntry -> Text)  -- ^ Custom format function
