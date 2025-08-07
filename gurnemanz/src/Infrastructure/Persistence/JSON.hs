{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE DefaultSignatures #-}

module Infrastructure.Persistence.JSON
    ( -- * JSON conversion
      ToJSON(..)
    , FromJSON(..)
    , ToJSONKey(..)
    , FromJSONKey(..)
    , ToCustomJSON(..)
      -- * Conversion helpers
    , encodeJSON
    , decodeJSON
    , encodeJSONFile
    , decodeJSONFile
    , encodeBinary
    , decodeBinary
    ) where

import Data.Aeson
    ( ToJSON(..), FromJSON(..)
    , ToJSONKey(..), FromJSONKey(..)
    , encode, decode, eitherDecode
    , Value(..)
    )
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy as BL
import qualified Data.ByteString as BS
import qualified Data.ByteString.Base64 as B64
import Data.Text (Text)
import qualified Data.Text as T
import Data.Text.Encoding (encodeUtf8, decodeUtf8)
import Control.Exception (try, SomeException)
import Control.Monad.IO.Class (MonadIO, liftIO)
import Data.Bifunctor (bimap)

import Core.Error (Error(..))
import Domain.Molecule.Types ()

-- | Type class for custom JSON encoding
class ToCustomJSON a where
    toCustomJSON :: a -> Value
    default toCustomJSON :: ToJSON a => a -> Value
    toCustomJSON = toJSON

-- | Helper for JSON encoding
encodeJSON :: ToJSON a => a -> ByteString
encodeJSON = encode

-- | Helper for JSON decoding
decodeJSON :: FromJSON a => ByteString -> Either String a
decodeJSON = eitherDecode

-- | Save JSON to file
encodeJSONFile :: (MonadIO m, ToJSON a) => FilePath -> a -> m (Either Error ())
encodeJSONFile path val = liftIO $ do
    result <- try $ BL.writeFile path (encode val)
    pure $ case result of
        Left e -> Left $ DataError $ "Failed to write JSON: " <> T.pack (show (e :: SomeException))
        Right _ -> Right ()

-- | Load JSON from file
decodeJSONFile :: (MonadIO m, FromJSON a) => FilePath -> m (Either Error a)
decodeJSONFile path = liftIO $ do
    result <- try $ BL.readFile path
    pure $ case result of
        Left e -> Left $ DataError $ "Failed to read JSON: " <> T.pack (show (e :: SomeException))
        Right contents ->
            case eitherDecode contents of
                Left e -> Left $ DataError $ "Failed to parse JSON: " <> T.pack e
                Right val -> Right val

-- | Handle binary data encoding
encodeBinary :: ByteString -> Text
encodeBinary = decodeUtf8 . B64.encode . BS.concat . BL.toChunks

-- | Handle binary data decoding
decodeBinary :: Text -> Either Text ByteString
decodeBinary t = BL.fromStrict <$> bimap T.pack id (B64.decode $ encodeUtf8 t)
