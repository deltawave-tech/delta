{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}

module Core.Lazy
    ( -- * Lazy computation wrapper
      Lazy(..)  -- Export constructor and fields
    , forceLazy
    , createLazy
    , fromLazy  -- Add export
    , toLazy    -- Add export
    ) where

import Control.Monad.IO.Class (MonadIO, liftIO)
import Data.IORef (IORef, newIORef, readIORef, writeIORef)
import Data.Aeson (ToJSON(..), FromJSON(..), object, Value)
import System.IO.Unsafe (unsafePerformIO)

-- | Lazy computation wrapper
data Lazy a = Lazy {
    _computed :: IORef (Maybe a),
    _computation :: IO a
}

-- | Show instance that doesn't force evaluation
instance Show (Lazy a) where
    show _ = "Lazy"

-- | JSON instances
instance ToJSON a => ToJSON (Lazy a) where
    toJSON = toJSON . fromLazy

instance FromJSON a => FromJSON (Lazy a) where
    parseJSON = fmap toLazy . parseJSON

-- | Force evaluation of a lazy computation
forceLazy :: MonadIO m => Lazy a -> m a
forceLazy (Lazy ref comp) = liftIO $ do
    mv <- readIORef ref
    case mv of
        Just v -> pure v
        Nothing -> do
            v <- comp
            writeIORef ref (Just v)
            pure v

-- | Create a new lazy computation
createLazy :: MonadIO m => IO a -> m (Lazy a)
createLazy comp = liftIO $ do
    ref <- newIORef Nothing
    pure $ Lazy {
        _computed = ref,
        _computation = comp
    }

-- | Extract the value from a Lazy computation
fromLazy :: Lazy a -> a
fromLazy l = unsafePerformIO $ forceLazy l

-- | Create a Lazy computation from a value
toLazy :: a -> Lazy a
toLazy a = Lazy {
    _computed = unsafePerformIO $ newIORef (Just a),
    _computation = pure a
}
