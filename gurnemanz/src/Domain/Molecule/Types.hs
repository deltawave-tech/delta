{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE OverloadedStrings #-}

module Domain.Molecule.Types
    ( -- * Core types
      ValidatedMolecule(..)
      -- * Re-exports
    , Valid
    , Draft
    -- * Smart constructors
    , createEmptyMolecule
    -- * Operations
    , storeMolecule
    ) where

import GHC.Generics (Generic)
import Data.Typeable (Typeable)
import Data.Text (Text)
import Data.Map (Map)
import qualified Data.Map as Map
import Data.Aeson (ToJSON(..), FromJSON(..), object, (.=), (.:), (.:?), (.!=), withObject, Value)
import Control.Monad.IO.Class (MonadIO)

import Core.Lazy (Lazy)
import Domain.Base.Types (Valid, Draft)

-- | A validated molecule with all its properties
data ValidatedMolecule a = ValidatedMolecule
    { moleculeId :: Text             -- Unique identifier for the molecule
    , parentMoleculeId :: Maybe Text -- NEW: The parent molecule this was derived from
    , transformationId :: Maybe Text -- NEW: The transformation that created this molecule
    , rationale :: Maybe Text        -- NEW: Captures the reason for the molecule creation/modification
    , validationInfo :: Maybe Text
    , smiles :: Maybe Text
    , inchiKey :: Maybe Text
    , molecularFormula :: Maybe Text
    , molecularWeight :: Maybe Double
    , xlogp :: Maybe Double          -- Maps to "logp"
    , hbd :: Maybe Int
    , hba :: Maybe Int
    , tpsa :: Maybe Double
    , rb :: Maybe Int
    , qed :: Maybe Double
    , sas :: Maybe Double            -- Maps to "sa_score" (assuming synonymy; adjust if distinct)
    , lipinskiCount :: Maybe Double
    , s3Path :: Maybe Text
    , dockingScore :: Maybe Double   -- Maps to "docking"
    , favorite :: Maybe Bool
    , plipInteractions :: Maybe Text -- From paper's "plip_interactions"
    , status :: Maybe Text           -- e.g., "modified", "new"
    , activity :: Maybe Text         -- NEW: Indicates if the molecule is "active" or "inactive" on a target
    , properties :: Map Text Value   -- Additional properties from external agents
    } deriving (Show, Generic, Typeable)

-- | Create an empty validated molecule
createEmptyMolecule :: Text -> ValidatedMolecule a
createEmptyMolecule id = ValidatedMolecule
    { moleculeId = id
    , parentMoleculeId = Nothing
    , transformationId = Nothing
    , rationale = Nothing
    , validationInfo = Nothing
    , smiles = Nothing
    , inchiKey = Nothing
    , molecularFormula = Nothing
    , molecularWeight = Nothing
    , xlogp = Nothing
    , hbd = Nothing
    , hba = Nothing
    , tpsa = Nothing
    , rb = Nothing
    , qed = Nothing
    , sas = Nothing
    , lipinskiCount = Nothing
    , s3Path = Nothing
    , dockingScore = Nothing
    , favorite = Nothing
    , plipInteractions = Nothing
    , status = Nothing
    , activity = Nothing
    , properties = Map.empty
    }

-- | Store a molecule (simplified replacement for what was in Operations)
storeMolecule :: MonadIO m => ValidatedMolecule Valid -> m ()
storeMolecule molecule = do
    -- In a real implementation, this would persist the molecule to storage
    -- For now, this is just a placeholder
    pure ()

-- | JSON instances
instance ToJSON (ValidatedMolecule a) where
    toJSON mol = object
        [ "moleculeId" .= moleculeId mol
        , "parentMoleculeId" .= parentMoleculeId mol
        , "transformationId" .= transformationId mol
        , "rationale" .= rationale mol
        , "validation" .= validationInfo mol
        , "smiles" .= smiles mol
        , "inchiKey" .= inchiKey mol
        , "molecularFormula" .= molecularFormula mol
        , "molecularWeight" .= molecularWeight mol
        , "xlogp" .= xlogp mol
        , "hbd" .= hbd mol
        , "hba" .= hba mol
        , "tpsa" .= tpsa mol
        , "rb" .= rb mol
        , "qed" .= qed mol
        , "sas" .= sas mol
        , "lipinskiCount" .= lipinskiCount mol
        , "s3Path" .= s3Path mol
        , "dockingScore" .= dockingScore mol
        , "favorite" .= favorite mol
        , "plipInteractions" .= plipInteractions mol
        , "status" .= status mol
        , "activity" .= activity mol
        , "properties" .= properties mol
        ]

instance FromJSON (ValidatedMolecule a) where
    parseJSON = withObject "ValidatedMolecule" $ \v ->
        ValidatedMolecule
            <$> v .: "moleculeId"
            <*> v .:? "parentMoleculeId"
            <*> v .:? "transformationId"
            <*> v .:? "rationale"
            <*> v .:? "validation"
            <*> v .:? "smiles"
            <*> v .:? "inchiKey"
            <*> v .:? "molecularFormula"
            <*> v .:? "molecularWeight"
            <*> v .:? "xlogp"
            <*> v .:? "hbd"
            <*> v .:? "hba"
            <*> v .:? "tpsa"
            <*> v .:? "rb"
            <*> v .:? "qed"
            <*> v .:? "sas"
            <*> v .:? "lipinskiCount"
            <*> v .:? "s3Path"
            <*> v .:? "dockingScore"
            <*> v .:? "favorite"
            <*> v .:? "plipInteractions"
            <*> v .:? "status"
            <*> v .:? "activity"
            <*> v .:? "properties" .!= Map.empty
