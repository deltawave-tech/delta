{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}

module Domain.Base.Graph
    ( -- * Re-export all types
      module Domain.Base.Graph.Types
      -- * Graph operations
    , validateGraph
    ) where

import Domain.Base.Graph.Types
import Core.Error (Error)

-- | Validate a graph structure
validateGraph :: MoleculeGraph -> Either Error MoleculeGraph
validateGraph = undefined  -- Implement validation logic
