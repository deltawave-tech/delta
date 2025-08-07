{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE OverloadedStrings #-}

module App.Types
    ( -- * Application Monad
      App(..)
    , runApp
      -- * Environment
    , AppEnv(..)
    ) where

import Control.Monad.Reader (ReaderT, MonadReader(..), runReaderT, lift)
import Control.Monad.Except (ExceptT, MonadError(..), runExceptT, throwError)
import Control.Monad.IO.Class (MonadIO, liftIO)
import qualified Data.UUID.V4 as UUID
import Data.UUID (UUID)
import qualified Data.Text as T

import Core.Error (Error)
import Core.Monad (MonadChem(..))
import Core.Base
    ( MonadLogger(..)
    , UserId
    , MoleculeId
    , SessionId
    )
import Common.Types
    ( Resources(..)
    , HasResources(..)
    , HasConfig(..)
    , LogHandle(..)
    , AppConfig(..) 
    , LogConfig(..)
    , LogContext(..)
    , LogLevel(..)
    )
import Infrastructure.Persistence.Logging
    ( LoggingT(..)
    , logMessage
    , runLoggingT
    )
import Domain.Molecule.Types (ValidatedMolecule, Valid, createEmptyMolecule)
import Domain.Base.Graph (MoleculeGraph)
import Domain.Transform.Types (Transformation, transformationId)
import qualified Domain.Transform.Operations as Trans

-- | Application environment
data AppEnv = AppEnv
    { appConfig :: !AppConfig
    , appResources :: !Resources
    }

-- | Application monad
newtype App a = App
    { unApp :: LoggingT (ReaderT AppEnv (ExceptT Error IO)) a }
    deriving newtype
        ( Functor
        , Applicative
        , Monad
        , MonadIO
        )

instance MonadError Error App where
    throwError e = App $ LoggingT $ \_ _ _ -> lift $ throwError e
    catchError (App (LoggingT m)) h = App $ LoggingT $ \handle config ctx ->
        catchError (m handle config ctx) (\e ->
            runLoggingT (unApp $ h e) handle config ctx)

instance MonadReader AppEnv App where
    ask = App $ LoggingT $ \_ _ _ -> ask
    local f (App (LoggingT m)) = App $ LoggingT $ \h c ctx ->
        local f (m h c ctx)

instance MonadLogger App where
    logDebug msg = App $ LoggingT $ \handle config ctx -> 
        liftIO $ putStrLn $ "[DEBUG] " ++ T.unpack msg

    logInfo msg = App $ LoggingT $ \handle config ctx -> 
        liftIO $ putStrLn $ "[INFO] " ++ T.unpack msg

    logWarning msg = App $ LoggingT $ \handle config ctx -> 
        liftIO $ putStrLn $ "[WARNING] " ++ T.unpack msg

    logError msg = App $ LoggingT $ \handle config ctx -> 
        liftIO $ putStrLn $ "[ERROR] " ++ T.unpack msg


instance MonadChem App where
    validateMolecule _ = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[CHEM] Validating molecule"
        pure $ createEmptyMolecule "placeholder-id"

    applyTransform _ _ = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[CHEM] Applying transformation"
        pure $ Right $ createEmptyMolecule "placeholder-id"

    rollbackTransformation tid chain = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[CHEM] Rolling back transformation"
        pure chain
            
    startSession uid _ = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[SESSION] Starting new session"
        uuid <- liftIO UUID.nextRandom
        pure uuid
        
    logTransformation trans = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[SESSION] Logging transformation"
        pure $ transformationId trans
        
    getTransformationChain molId = App $ LoggingT $ \handle config ctx -> do
        liftIO $ putStrLn "[SESSION] Retrieving transformation chain"
        pure []

-- | Environment instances
instance HasConfig AppEnv where
    getConfig = appConfig

instance HasResources AppEnv where
    getResources = appResources

-- | Run the application monad
runApp :: AppEnv -> App a -> IO (Either Error a)
runApp env (App action) =
    runExceptT $
        flip runReaderT env $
            runLoggingT action
                (resourceLogHandle $ appResources env)
                (logConfig $ appConfig env)
                defaultLogContext

-- | Default logging context
defaultLogContext :: LogContext
defaultLogContext = LogContext
    { contextLevel = Info
    , contextComponent = "App"
    , contextCorrelationId = Nothing
    , contextTags = []
    }