{-# LANGUAGE GADTs #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE FlexibleContexts #-}

module Domain.Session.Types
    ( -- * Core types
      SessionF(..)
      -- * Smart constructors
    , startSession
    , logTransformation
    , getTransformationChain
    , rollbackTo
    ) where

import GHC.Generics (Generic)
import Control.Monad.IO.Class (MonadIO)
import Control.Monad.Error.Class (MonadError)
import Data.Text (Text)

import Core.Base (UserId, SessionId, MolId, MoleculeId)
import Core.Error (Error)
import Domain.Transform.Types (TransformationId, TransformationChain, Transformation)

-- | Session free monad for operations
data SessionF next where
    StartSession :: UserId -> Transformation -> (SessionId -> next) -> SessionF next
    LogTransformation :: Transformation -> (TransformationId -> next) -> SessionF next
    GetTransformationChain :: MoleculeId -> ([Transformation] -> next) -> SessionF next
    RollbackTo :: TransformationId -> (TransformationChain -> next) -> SessionF next

-- Smart constructors
startSession :: MonadIO m => UserId -> Transformation -> m SessionId
startSession = undefined

logTransformation :: MonadIO m => Transformation -> m TransformationId
logTransformation = undefined

getTransformationChain :: MonadIO m => MoleculeId -> m [Transformation]
getTransformationChain = undefined

rollbackTo :: (MonadIO m, MonadError Error m) => TransformationId -> m TransformationChain
rollbackTo = undefined
