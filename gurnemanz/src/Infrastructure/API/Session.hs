{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}

module Infrastructure.API.Session
    ( handleSessionRequest
    , SessionRequest(..)
    , SessionResponse(..)
    ) where

import qualified Data.UUID.V4 as UUID
import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Error.Class (MonadError, throwError)
import Data.Text (Text)
import qualified Data.Text as T
import Data.List.NonEmpty (NonEmpty(..))

import Core.Error (Error(..))
import Core.Base
    ( MonadLogger(..)
    , SessionId
    , UserId
    , TransformationId
    )
import Domain.Transform.Types
    ( TransformationChain(..)
    , Transformation(..)
    , createTransformation
    , createEmptyChain
    )
import Domain.Molecule.Types (ValidatedMolecule, Valid)
import qualified Domain.Transform.Operations as Trans
import Core.Monad (MonadChem(..))

-- | Handle session operations
handleSessionRequest
    :: (MonadIO m, MonadError Error m, MonadLogger m, MonadChem m)
    => SessionRequest
    -> m SessionResponse
handleSessionRequest = \case
    CreateNewSession _ -> do
        sessionId <- liftIO UUID.nextRandom
        initTrans <- createTransformation "session-init"
        let systemTrans = initTrans { agent = "system" }
        let _ = createEmptyChain (systemTrans :| [])

        logInfo $ "Created new session: " <> T.pack (show sessionId)
        pure $ SessionCreated sessionId

    GetExistingSession _ -> do
        throwError $ ValidationError "Session not found"

    ApplyTransformToSession _ _ _ -> do
        throwError $ ValidationError "Session not found"

    RollbackSession _ _ -> do
        throwError $ ValidationError "Session not found"

-- | Session request types
data SessionRequest
    = CreateNewSession UserId
    | GetExistingSession SessionId
    | ApplyTransformToSession SessionId Text (ValidatedMolecule Valid)
    | RollbackSession SessionId TransformationId

-- | Session response types
data SessionResponse
    = SessionCreated SessionId
    | SessionRetrieved TransformationChain
    | TransformationApplied (ValidatedMolecule Valid) TransformationChain
    | SessionRolledBack TransformationChain