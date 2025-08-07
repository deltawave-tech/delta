{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

module Domain.Base.Property
    ( -- * Property types
      Property(..)
    , PropertyId(..)
      -- * JSON instances
    , module Data.Aeson
    ) where

import GHC.Generics (Generic)
import Data.Text (Text)
import Data.Aeson
    ( ToJSON(..), FromJSON(..)
    , ToJSONKey(..), FromJSONKey(..)
    , Value(String), object, withText
    )
import Data.Aeson.Types (ToJSONKeyFunction(..), FromJSONKeyFunction(..))
import Data.Aeson.Encoding (text)

-- | Property value container
data Property = Property
    deriving (Show, Generic)

-- | Property identifiers
data PropertyId
    = LogPProperty
    | TopologyProperty
    deriving (Eq, Ord, Show, Generic)

-- | JSON instances
instance ToJSON Property where
    toJSON _ = object []

instance FromJSON Property where
    parseJSON _ = pure Property

instance ToJSON PropertyId where
    toJSON = \case
        LogPProperty -> String "LogP"
        TopologyProperty -> String "Topology"

instance FromJSON PropertyId where
    parseJSON = withText "PropertyId" $ \case
        "LogP" -> pure LogPProperty
        "Topology" -> pure TopologyProperty
        _ -> fail "Invalid PropertyId"

instance ToJSONKey PropertyId where
    toJSONKey = ToJSONKeyText propertyIdToText propertyIdToEncodingText
      where
        propertyIdToText LogPProperty = "LogP"
        propertyIdToText TopologyProperty = "Topology"
        propertyIdToEncodingText LogPProperty = text "LogP"
        propertyIdToEncodingText TopologyProperty = text "Topology"

instance FromJSONKey PropertyId where
    fromJSONKey = FromJSONKeyText $ \case
        "LogP" -> LogPProperty
        "Topology" -> TopologyProperty
        _ -> error "Invalid PropertyId"
