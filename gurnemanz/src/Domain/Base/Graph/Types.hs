module Domain.Base.Graph.Types
    ( -- * Basic types
      Atom(..)
    , Bond(..)
    , Graph(..)
    , BaseGraph
      -- * Properties
    , Properties(..)
      -- * Enhanced graph
    , MoleculeGraph(..)
    , createEmptyGraph
      -- * Initialization
    , initialGraph
    , initialProps
    ) where

import Data.Map (Map)
import qualified Data.Map as Map
import GHC.Generics (Generic)
import Data.Aeson (ToJSON(..), FromJSON(..), (.=), (.:), (.:?), withObject, object)
import Data.Aeson.Types ((.!=))

import Core.Lazy (Lazy)
import Domain.Base.Property (PropertyId, Property)

-- Basic types
data Atom = Atom
    deriving stock (Show, Generic)
    deriving anyclass (ToJSON, FromJSON)

data Bond = Bond
    deriving stock (Show, Generic)
    deriving anyclass (ToJSON, FromJSON)

data Graph v e = Graph
    deriving stock (Show, Generic)

instance FromJSON (Graph v e) where
    parseJSON = withObject "Graph" $ \_ -> pure Graph

instance ToJSON (Graph v e) where
    toJSON _ = object []

-- | Properties container
data Properties = Properties
    deriving stock (Show, Generic)

instance FromJSON Properties where
    parseJSON = withObject "Properties" $ \_ -> pure Properties

instance ToJSON Properties where
    toJSON _ = object []

type BaseGraph = Graph Atom Bond

-- Enhanced graph
data MoleculeGraph = MoleculeGraph
    { _structure :: BaseGraph
    , _properties :: Properties
    , _computedProperties :: Map PropertyId (Lazy Property)
    } deriving stock (Show, Generic)

instance FromJSON MoleculeGraph where
    parseJSON = withObject "MoleculeGraph" $ \v -> do
        struct <- v .: "structure"
        props <- v .: "properties"
        comps <- v .:? "computedProperties" .!= Map.empty
        pure $ MoleculeGraph
            { _structure = struct
            , _properties = props
            , _computedProperties = comps
            }

instance ToJSON MoleculeGraph where
    toJSON mg = object
        [ "structure" .= _structure mg
        , "properties" .= _properties mg
        , "computedProperties" .= _computedProperties mg
        ]

createEmptyGraph :: MoleculeGraph
createEmptyGraph = MoleculeGraph
    { _structure = Graph
    , _properties = Properties
    , _computedProperties = Map.empty
    }

-- | Create initial empty graph
initialGraph :: BaseGraph
initialGraph = Graph

-- | Create initial empty properties
initialProps :: Properties
initialProps = Properties
