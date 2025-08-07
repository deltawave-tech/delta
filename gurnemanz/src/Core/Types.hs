-- src/Core/Types.hs
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE CPP #-}

module Core.Types
    ( -- * Resource Types (re-exported)
      Resources(..)
      -- * Core Classes (re-exported)
    , HasResources(..)
    , HasConfig(..)
    -- * Other types (re-exported)
    , LogHandle(..)
    , AppConfig(..)
    ) where

-- Import everything from Common.Types
import Common.Types
    ( Resources(..)
    , HasResources(..)
    , HasConfig(..)
    , LogHandle(..)
    , AppConfig(..)
    )