{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}

module Domain.Transform.Operations
    ( -- * Chain operations
      rollbackTo
    , getChainForMolecule
    , addTransformation
    , addMolecules  -- New function to add molecules to a chain
      -- * Validation
    , validateTransformation
      -- * Extended operations
    , createInteraction
    , updateTransformation
    ) where

import Control.Monad.IO.Class (MonadIO)
import Control.Monad.Error.Class (MonadError, throwError)
import Data.Map (Map)
import qualified Data.Map as Map
import qualified Data.List.NonEmpty as NE
import Data.Text (Text)
import qualified Data.Text as Text
import Data.Aeson (Value)
import qualified Data.Aeson as Aeson

import Core.Error (Error(..))
import Core.Base (Result, MonadLogger(..), MoleculeId)
import Core.Lazy (Lazy)
import Domain.Transform.Types
import Domain.Molecule.Types (ValidatedMolecule, Valid)

-- | Create a new interaction transformation
createInteraction
    :: (MonadIO m, MonadError Error m, MonadLogger m)
    => Text         -- ^ User message
    -> Maybe MoleculeId  -- ^ Input molecule
    -> m Transformation
createInteraction msg molId = do
    logInfo $ "Creating interaction with message: " <> msg
    trans <- createTransformation "user-interaction"
    pure $ trans
        { userMessage = Just msg
        , inputMoleculeIds = maybe [] (:[]) molId
        , agent = "user" -- Set a default agent
        }

-- | Update transformation with agent response
updateTransformation
    :: (MonadError Error m, MonadLogger m)
    => Transformation
    -> Text  -- ^ Agent response
    -> Maybe MethodDetails
    -> [MoleculeId]  -- ^ Output molecules
    -> m Transformation
updateTransformation trans response method outputs = do
    logInfo $ "Updating transformation: " <> Text.pack (show $ transformationId trans)
    pure $ trans
        { agentResponse = Just response
        , methodDetails = method
        , outputMoleculeIds = outputs
        }

-- | Roll back a transformation chain to a specific point
rollbackTo
    :: (MonadError Error m, MonadLogger m)
    => TransformationId
    -> TransformationChain
    -> m TransformationChain
rollbackTo tid chain = do
    logInfo $ "Rolling back transformation: " <> Text.pack (show tid)
    let steps = NE.toList (_steps chain)
        (_, rest) = break (\s -> transformationId s == tid) steps
    case rest of
        [] -> throwError $ TransformationError "Target transformation not found"
        (target:_) -> case NE.nonEmpty [target] of
            Nothing -> throwError $ TransformationError "Failed to create valid chain"
            Just validNE -> pure $ TransformationChain
                { _steps = validNE
                , _status = RolledBack
                , _molecules = _molecules chain  -- Preserve molecule data
                }


-- | Validate a transformation
validateTransformation
    :: (MonadError Error m)
    => Transformation
    -> m ()
validateTransformation trans =
    if Text.null (transformationType trans)
        then throwError $ ValidationError "Empty transformation type"
        else pure ()

-- | Get the chain of transformations leading to a specific molecule
getChainForMolecule :: TransformationChain -> MoleculeId -> [Transformation]
getChainForMolecule chain molId = reverse $ go molId []
  where
    steps = NE.toList (_steps chain)
    go currentId acc =
      case filter (\t -> currentId `elem` outputMoleculeIds t) steps of
        [] -> acc
        (t:_) -> case inputMoleculeIds t of
            [] -> t:acc
            (parentId:_) -> go parentId (t:acc)

-- | Add molecules to a transformation chain
addMolecules
    :: TransformationChain
    -> [(MoleculeId, Value)]  -- ^ List of molecule IDs and their data
    -> TransformationChain
addMolecules chain moleculePairs =
    chain { _molecules = Map.union (Map.fromList moleculePairs) (_molecules chain) }

-- Internal helpers
addTransformation
    :: TransformationChain
    -> Transformation
    -> TransformationChain
addTransformation chain trans =
    chain { _steps = NE.cons trans (_steps chain) }

-- Helper function
whenJust :: Applicative f => Maybe a -> (a -> f ()) -> f ()
whenJust Nothing _ = pure ()
whenJust (Just x) f = f x