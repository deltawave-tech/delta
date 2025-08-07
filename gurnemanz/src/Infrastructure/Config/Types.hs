{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE LambdaCase #-}

module Infrastructure.Config.Types
    ( -- * Core Configuration Types (re-exported)
      AppConfig(..)
    , Environment(..)
      -- * Component Configurations (re-exported)
    , ChemConfig(..)
    , APIConfig(..)
      -- * Type Re-exports
    , Word16
    ) where

-- Import from Common.Types
import Common.Types
    ( AppConfig(..)
    , Environment(..)
    , ChemConfig(..)
    , APIConfig(..)
    , Word16
    )