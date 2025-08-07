{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}

module Domain.Session.Interpreter
    ( -- * Session State
      SessionState(..)
      -- * Interpreter
    , interpretSessionF
    ) where

import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Error.Class (MonadError, throwError)
import Control.Monad.State (StateT, get, put)
import Control.Monad.Trans (lift)
import Data.List.NonEmpty (NonEmpty)
import qualified Data.List.NonEmpty as NE
import Data.UUID.V4 as UUID (nextRandom)
import Data.Text (Text)
import qualified Data.Text as Text

import Core.Base (UserId, SessionId, MoleculeId, TransformationId, MonadLogger(..))
import Core.Error (Error(..))
import Domain.Session.Types (SessionF(..))
import Domain.Transform.Types (TransformationChain, Transformation, transformationId, createEmptyChain)
import Domain.Transform.Operations (getChainForMolecule, addTransformation, rollbackTo)

-- | Session state
data SessionState = SessionState
    { sessionId :: SessionId
    , userId :: UserId
    , transformationChain :: TransformationChain
    }

-- | Interpreter for SessionF
interpretSessionF :: (MonadIO m, MonadError Error m, MonadLogger m) => SessionF a -> StateT SessionState m a
interpretSessionF (StartSession uid initialTrans next) = do
    sid <- liftIO UUID.nextRandom
    -- Create the chain with the initial transformation
    let initialChain = createEmptyChain (initialTrans NE.:| [])
    put SessionState
       { sessionId = sid
       , userId = uid
       , transformationChain = initialChain
       }
    pure $ next sid

interpretSessionF (LogTransformation trans next) = do
    state <- get
    let chain = transformationChain state
        newChain = addTransformation chain trans
    put state { transformationChain = newChain }
    pure $ next (transformationId trans)

interpretSessionF (GetTransformationChain molId next) = do
    state <- get
    let chain = transformationChain state
        transformations = getChainForMolecule chain molId
    pure $ next transformations

interpretSessionF (RollbackTo tid next) = do
    state <- get
    let chain = transformationChain state
    -- Use lift to apply the MonadLogger constraint from the base monad
    newChain <- lift $ rollbackTo tid chain 
    put state { transformationChain = newChain }
    pure $ next newChain