{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE CPP #-}

module Core.Error
    ( -- * Error Types
      Error(..)
    , AsError(..)
      -- * Error handling
    , throwError
    , catchError
    , handleError
    , mapError
    , runError
    , showError
    ) where

import GHC.Generics (Generic)
import Data.Text (Text)
import qualified Data.Text as Text
import Control.Monad.Except (MonadError(..), ExceptT, runExceptT)
import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Reader (ReaderT)
import Data.Aeson
    ( ToJSON(..), FromJSON(..), object, (.=), (.:)
    , withObject, fromJSON
    )
import Data.Aeson.Types (Result(..))
import Control.Exception (try, Exception)

import Core.Base (MonadLogger(..))


-- | Core error types
data Error
    = ValidationError Text
    | ComputationError Text
    | TransformationError Text
    | DataError Text
    deriving (Show, Eq, Generic, Exception)

-- | Class for converting to Error
class AsError e where
    -- | Convert a specific error type to our general Error type
    toError :: e -> Error
    -- | Convert our general Error type to a maybe of a specific type
    fromError :: Error -> Maybe e

-- | Instance for Error itself
instance AsError Error where
    toError = id
    fromError = Just


-- | JSON instances for Error
instance ToJSON Error where
    toJSON = \case
        ValidationError msg -> object ["type" .= ("validation" :: Text), "message" .= msg]
        ComputationError msg -> object ["type" .= ("computation" :: Text), "message" .= msg]
        TransformationError msg -> object ["type" .= ("transformation" :: Text), "message" .= msg]
        DataError msg -> object ["type" .= ("data" :: Text), "message" .= msg]

instance FromJSON Error where
    parseJSON = withObject "Error" $ \v -> do
        typ <- v .: "type"
        msg <- v .: "message"
        pure $ case (typ :: Text) of
            "validation" -> ValidationError msg
            "computation" -> ComputationError msg
            "transformation" -> TransformationError msg
            "data" -> DataError msg
            _ -> DataError "Unknown error type"


-- | Handle an error by converting it and logging
handleError :: (MonadLogger m, MonadError Error m, AsError e) => e -> m a
handleError err = do
    let errConverted = toError err
    logError $ errorMessage errConverted
    throwError errConverted
  where
    errorMessage :: Error -> Text
    errorMessage = showError

-- | Show error as text
showError :: Error -> Text
showError = \case
    ValidationError msg -> "Validation error: " <> msg
    ComputationError msg -> "Computation error: " <> msg
    TransformationError msg -> "Transformation error: " <> msg
    DataError msg -> "Data error: " <> msg

-- | Map one error type to another
mapError :: (AsError e1, AsError e2) => (e1 -> e2) -> Error -> Error
mapError f err = case fromError err of
    Just e1 -> toError (f e1)
    Nothing -> err

-- | Run an error handler in IO
runError :: MonadIO m => ExceptT Error m a -> m (Either Error a)
runError = runExceptT

