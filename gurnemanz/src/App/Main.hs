module Main where

import Control.Monad.IO.Class (liftIO)
import Control.Monad (forever)
import qualified Data.Text.IO as TIO
import Data.Text (Text)
import qualified Data.Text as T
import System.Exit (exitFailure)
import Control.Monad.Except (runExceptT)

import App.Types (App, runApp)
import App.Env (withAppEnv)
import Infrastructure.Config.Loader (loadConfig, validateConfig)
import Infrastructure.Config.Types (AppConfig(..))
import Infrastructure.Config.Environment (validateEnvironment, getEnvironmentConfig)
import Core.Error (Error(..), showError)
import Core.Base (MonadLogger(..))
import Domain.Molecule.Types (createEmptyMolecule)
import Domain.Transform.Types (createTransformation, Transformation(..))
import Domain.Session.Types (startSession, logTransformation)
import Infrastructure.Persistence.Logging.Console (withConsoleLogger)
import Infrastructure.API.Socket (runSocketServer)

-- | Main application logic
mainApp :: App ()
mainApp = do
    logInfo "Starting main application logic"
    
    -- Create an initial transformation for the session
    initialTrans <- liftIO $ createTransformation "session_start"
    
    -- Update the transformation fields manually
    let initialTransWithDetails = initialTrans 
          { agent = "system"
          }
    
    -- Start a session with the initial transformation
    sid <- startSession "test-user-id" initialTransWithDetails
    logInfo $ "Created session with ID: " <> showText sid
    
    -- Create a sample molecule directly
    let validated = createEmptyMolecule "test-molecule"
    
    
    -- Create a test transformation for molecule modification
    transformation <- liftIO $ createTransformation "molecule_modify"
    let transformationWithDetails = transformation
          { agent = "agent-1"
          }
    
    -- Log the transformation
    tid <- logTransformation transformationWithDetails
    logInfo $ "Logged transformation with ID: " <> showText tid
    
    logInfo "Main application logic completed"

-- | Show a value as Text
showText :: Show a => a -> Text
showText = T.pack . show

-- | Error handling for the main application
handleAppError :: Error -> IO ()
handleAppError err = do
    TIO.putStrLn $ "Application error: " <> showError err
    exitFailure

-- | Main entry point
main :: IO ()
main = withConsoleLogger $ do
    putStrLn "GURNEMANZ Starting..."

    configResult <- runExceptT $ do
        baseConfig <- loadConfig "config.yaml"
        validConfig <- validateConfig baseConfig
        validateEnvironment validConfig
        let envConfig = getEnvironmentConfig (environment validConfig)
        pure $ envConfig { environment = environment validConfig }

    case configResult of
        Left err -> handleAppError $ ValidationError $ "Configuration error: " <> showError err
        Right config -> withAppEnv config $ \env -> do
            logInfo "Starting application with configuration..."
            -- Start socket server in the background
            runSocketServer env (apiConfig config)
            -- Run the main application
            result <- runApp env mainApp
            case result of
                Left err -> handleAppError err
                Right _ -> do
                    logInfo "Operation completed successfully"
                    forever $ pure ()  -- Keep the server running
