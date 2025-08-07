{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE CPP #-}

module Core.Monad
    ( -- * Type Classes
      MonadChem(..)
      -- * Actions
    , ChemAction(..)
    ) where

import Control.Monad.Reader (ReaderT(..), MonadReader, asks, lift)
import Control.Monad.Except (ExceptT, MonadError, runExceptT, throwError, catchError)
import Control.Monad.IO.Class (MonadIO, liftIO)

import Core.Error (Error(..))
import Core.Base
    ( MoleculeId
    , TransformationId
    , SessionId
    , UserId
    , MonadLogger(..)
    , LogLevel(..)
    )
import Common.Types
    ( LogConfig(..)
    , LogContext(..)
    , LogHandle(..)
    )

import Domain.Molecule.Types (ValidatedMolecule, Valid)
import Domain.Transform.Types (Transformation, TransformationChain)
import Infrastructure.Persistence.Logging (LoggingT(..), logMessage)
import qualified Data.Aeson as Aeson
import Data.Aeson (Value)
import Data.Text (Text)

-- | Session and chemistry capabilities
class (MonadError Error m) => MonadChem m where
    validateMolecule :: String -> m (ValidatedMolecule Valid)  -- Changed from MoleculeGraph
    applyTransform :: Transformation -> ValidatedMolecule Valid -> m (Either Error (ValidatedMolecule Valid))
    rollbackTransformation :: TransformationId -> TransformationChain -> m TransformationChain
    startSession :: UserId -> Transformation -> m SessionId
    logTransformation :: Transformation -> m TransformationId
    getTransformationChain :: MoleculeId -> m [Transformation]


-- | Chemistry action monad
newtype ChemAction m a = ChemAction
    { runChem :: ExceptT Error m a }
    deriving newtype (Functor, Applicative, Monad, MonadIO)

