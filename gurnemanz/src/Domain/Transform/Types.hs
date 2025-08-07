{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}

module Domain.Transform.Types
    ( -- * Core types
      Transformation(..)
    , TransformationChain(..)
    , ChainStatus(..)
    , TransformationId
    , MethodDetails(..)
      -- * Chain constructors
    , createEmptyChain
    -- * Smart constructors
    , createTransformation
    -- * Re-exports
    , NonEmpty(..)
    ) where

import GHC.Generics (Generic)
import Data.Text (Text)
import Data.UUID (UUID)
import qualified Data.UUID.V4 as UUID
import Data.Map (Map)
import qualified Data.Map as Map
import Data.List.NonEmpty (NonEmpty(..))
import Control.Monad.IO.Class (MonadIO, liftIO)
import Data.Aeson (ToJSON(..), FromJSON(..), object, (.=), (.:), (.:?), (.!=), withObject, Value(..))

import Core.Base (TransformationId, Result, MoleculeId, MethodId)
import Core.Lazy (Lazy)

-- | Method application details
data MethodDetails = MethodDetails
    { methodId :: MethodId
    , parameters :: Map Text Value
    } deriving (Show, Generic)

-- | Status of a transformation chain
data ChainStatus
    = RolledBack
    | Active
    deriving (Show, Generic)

-- | A single transformation step
data Transformation = Transformation
    { transformationId :: TransformationId
    , transformationType :: Text
    , userMessage :: Maybe Text
    , agentResponse :: Maybe Text
    , methodDetails :: Maybe MethodDetails
    , inputMoleculeIds :: [MoleculeId]  -- Changed from Maybe MoleculeId to support multiple inputs
    , outputMoleculeIds :: [MoleculeId]
    , agent :: Text          -- Identifies the agent (e.g., "Medicinal Chemist", "AI_expert")
    , proteinSequence :: Maybe Text    -- Stores a protein sequence associated with the transformation
    , proteinPath :: Maybe Text        -- Stores a file path to the protein structure file
    , iteration :: Maybe Int           -- Stores the iteration number for this transformation
    } deriving (Show, Generic)

-- | A chain of transformations and their associated molecules
data TransformationChain = TransformationChain
    { _steps :: NonEmpty Transformation
    , _status :: ChainStatus
    , _molecules :: Map.Map MoleculeId Value  -- Store the original molecule data by ID
    } deriving (Show, Generic)

-- | Create a new transformation with simplified interface
createTransformation :: MonadIO m => Text -> m Transformation
createTransformation transformType = do
    tid <- liftIO UUID.nextRandom
    pure $ Transformation
        { transformationId = tid
        , transformationType = transformType
        , userMessage = Nothing
        , agentResponse = Nothing
        , methodDetails = Nothing
        , inputMoleculeIds = []
        , outputMoleculeIds = []
        , agent = "system"
        , proteinSequence = Nothing
        , proteinPath = Nothing
        , iteration = Nothing
        }

-- | Create an empty transformation chain
createEmptyChain :: NonEmpty Transformation -> TransformationChain
createEmptyChain steps = TransformationChain
    { _steps = steps
    , _status = Active
    , _molecules = Map.empty
    }

-- JSON instances
instance ToJSON ChainStatus
instance FromJSON ChainStatus

instance ToJSON MethodDetails where
    toJSON md = object
        [ "methodId" .= methodId md
        , "parameters" .= parameters md
        ]

instance FromJSON MethodDetails where
    parseJSON = withObject "MethodDetails" $ \v ->
        MethodDetails
            <$> v .: "methodId"
            <*> v .: "parameters"

instance ToJSON Transformation where
    toJSON t = object
        [ "transformationId" .= transformationId t
        , "transformationType" .= transformationType t
        , "userMessage" .= userMessage t
        , "agentResponse" .= agentResponse t
        , "methodDetails" .= methodDetails t
        , "inputMoleculeIds" .= inputMoleculeIds t
        , "outputMoleculeIds" .= outputMoleculeIds t
        , "agent" .= agent t
        , "proteinSequence" .= proteinSequence t
        , "proteinPath" .= proteinPath t
        , "iteration" .= iteration t
        ]

instance FromJSON Transformation where
    parseJSON = withObject "Transformation" $ \v ->
        Transformation
            <$> v .: "transformationId"
            <*> v .: "transformationType"
            <*> v .:? "userMessage"
            <*> v .:? "agentResponse"
            <*> v .:? "methodDetails"
            <*> v .:? "inputMoleculeIds" .!= []
            <*> v .: "outputMoleculeIds"
            <*> v .: "agent"
            <*> v .:? "proteinSequence"
            <*> v .:? "proteinPath"
            <*> v .:? "iteration"

instance ToJSON TransformationChain where
    toJSON chain = object
        [ "steps" .= _steps chain
        , "status" .= _status chain
        , "results" .= object [ "molecules" .= _molecules chain ]  -- Include molecules in results
        ]

instance FromJSON TransformationChain where
    parseJSON = withObject "TransformationChain" $ \v -> do
        steps <- v .: "steps"
        status <- v .: "status"
        -- Try to get molecules from results.molecules, default to empty map
        results <- v .:? "results" .!= object []
        molecules <- case results of
            Object o -> o .:? "molecules" .!= Map.empty
            _ -> pure Map.empty
        pure $ TransformationChain steps status molecules