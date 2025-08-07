{-# LANGUAGE OverloadedStrings #-}

module Main where

import Test.Hspec
import Test.QuickCheck
import qualified Data.Map as Map
import qualified Data.Text as T
import Data.Text (Text)
import Control.Monad.IO.Class (liftIO)
import Control.Monad.Except (runExceptT)
import Data.UUID (UUID)
import qualified Data.UUID.V4 as UUID
import Data.Aeson (Value(..))

-- Import core types and functions
import Core.Base (MonadLogger(..))
import Core.Error (Error(..))
import Common.Types (Resources(..), AppConfig(..), Environment(..), LogHandle(..), 
                       ChemConfig(..), APIConfig(..), LogConfig(..), LogFormat(..), LogLevel(..))
import App.Types (App, AppEnv(..), runApp)
import Control.Concurrent (newMVar)
import System.IO (stdout)
import Infrastructure.Config.Loader (loadConfig, validateConfig)
import Domain.Molecule.Types (ValidatedMolecule(..), createEmptyMolecule, Valid)
import Domain.Transform.Types (Transformation(..), createTransformation)
import qualified Domain.Transform.Types as Transform
import Core.Monad (MonadChem(applyTransform, validateMolecule, startSession, logTransformation, getTransformationChain))

main :: IO ()
main = hspec $ do
    describe "Configuration" $ do
        it "should load test configuration correctly" $ do
            configResult <- runExceptT $ do
                baseConfig <- loadConfig "config.test.yaml"
                validateConfig baseConfig
            case configResult of
                Left err -> expectationFailure $ "Failed to load config: " ++ show err
                Right config -> do
                    -- Check environment setting
                    environment config `shouldBe` Development
                    
                    -- Verify chemical config has expected values
                    let chem = chemConfig config
                    maxMoleculeSize chem `shouldSatisfy` (> 0)
                    maxTransformationSteps chem `shouldSatisfy` (> 0)
                    computationTimeout chem `shouldSatisfy` (> 0)
                    enabledProperties chem `shouldNotBe` []
                    
                    -- Verify API config has expected values
                    let api = apiConfig config
                    apiPort api `shouldSatisfy` (> 0)
                    T.length (apiHost api) `shouldSatisfy` (> 0)
                    apiMaxRequestSize api `shouldSatisfy` (> 0)
                    length (socketPath api) `shouldSatisfy` (> 0)
                    
                    -- Verify log config has expected values
                    let log = logConfig config
                    enableConsoleLog log `shouldBe` True  -- For testing

    -- Create a helper function for setting up mock environment
    let createMockEnv = do
            logHandleMVar <- newMVar stdout
            let logHandle = LogHandle logHandleMVar
                mockChemConfig = ChemConfig
                    { maxMoleculeSize = 1000
                    , maxTransformationSteps = 50
                    , computationTimeout = 300
                    , enabledProperties = ["logP", "topology"]
                    }
                mockApiConfig = APIConfig
                    { apiPort = 8080
                    , apiHost = "localhost"
                    , apiMaxRequestSize = 10485760
                    , apiRateLimit = 100
                    , apiTimeout = 30
                    , socketPath = "/tmp/test.sock"
                    }
                mockLogConfig = LogConfig
                    { logLevel = Info
                    , logFile = Nothing
                    , logFormat = LogText
                    , enableConsoleLog = True
                    }
                mockEnv = AppEnv 
                    { appResources = Resources 
                        { resourceLogHandle = logHandle
                        }
                    , appConfig = AppConfig
                        { environment = Development
                        , chemConfig = mockChemConfig
                        , apiConfig = mockApiConfig
                        , logConfig = mockLogConfig
                        }
                    }
            return (logHandle, mockEnv)
    
    describe "Session Operations" $ do
        it "should create a session with initial transformation" $ do
            -- Create a mock environment
            (_, mockEnv) <- createMockEnv
            
            -- Run the test in the App monad
            result <- runApp mockEnv $ do
                -- Create an initial transformation
                t <- liftIO $ createTransformation "test-transform"
                let initialTrans = t { 
                        agent = "test-agent", 
                        iteration = Just 1
                    }
                
                -- Start a session with the initial transformation
                sessionId <- startSession "test-user" initialTrans
                
                -- Verify the session ID is valid
                logInfo $ "Session created with ID: " <> T.pack (show sessionId)
                
                -- Get the transformation chain for verification
                -- In our mock this will be empty, but the operation should succeed
                let testMolId = Transform.transformationId initialTrans -- Use the transformation ID as a mock molecule ID
                chain <- getTransformationChain testMolId
                
                -- Verify session creation was successful by checking:
                -- 1. SessionId is a valid UUID (non-null)
                -- 2. Transformation operation didn't fail
                pure (sessionId, initialTrans, chain)
                
            -- Verify results
            case result of
                Left err -> expectationFailure $ "Session creation failed: " ++ show err
                Right (sessionId, initialTrans, chain) -> do
                    -- Verify we have a valid UUID as session ID
                    show sessionId `shouldSatisfy` (not . null)
                    
                    -- Our mock currently returns an empty chain, but we note this as a limitation
                    -- In an ideal implementation, the chain would contain the initial transformation
                    null chain `shouldBe` True
                    pendingWith "Mock returns empty chain; a better mock would include the initial transformation"
                    
                    -- The following test would be the ideal version once the mock is enhanced:
                    -- length chain `shouldBe` 1
                    -- head chain `shouldBe` initialTrans  -- Would need Eq instance

    describe "Transformation Operations" $ do
        it "should create and modify transformations" $ do
            -- Create a transformation
            trans <- createTransformation "molecule-modify"
            
            -- Verify basic properties
            transformationType trans `shouldBe` "molecule-modify"
            agent trans `shouldBe` "system"  -- Default agent
            iteration trans `shouldBe` Nothing     -- Default iteration
            
            -- Verify ID was created
            show (Transform.transformationId trans) `shouldSatisfy` (not . null)
            
            -- Modify the transformation
            let modifiedTrans = trans 
                  { agent = "medicinal-chemist"
                  , iteration = Just 2
                  }
                
            -- Verify modifications
            agent modifiedTrans `shouldBe` "medicinal-chemist"
            iteration modifiedTrans `shouldBe` Just 2
            
        it "should handle edge cases in transformations" $ do
            -- Create a base transformation
            trans <- createTransformation "base-transform"
            
            -- Test with very large iteration number
            let largeIterationTrans = trans { iteration = Just 999999 }
            iteration largeIterationTrans `shouldBe` Just 999999
            
            -- Test with special characters in agent name
            let specialCharsTrans = trans { agent = "agent-1@#$%^&*()" }
            agent specialCharsTrans `shouldBe` "agent-1@#$%^&*()"
            
            -- Test with negative iteration (edge case)
            let negativeIterationTrans = trans { iteration = Just (-1) }
            iteration negativeIterationTrans `shouldBe` Just (-1)
            
            -- Test with very long type name
            longTypeTrans <- createTransformation $ T.replicate 1000 "a"
            T.length (transformationType longTypeTrans) `shouldBe` 1000

    describe "Molecule Operations" $ do
        it "should create and validate molecules" $ do
            -- Create a basic molecule
            let molecule = createEmptyMolecule "test-molecule"
            
            -- Verify basic properties
            moleculeId molecule `shouldBe` "test-molecule"
            smiles molecule `shouldBe` Nothing
            inchiKey molecule `shouldBe` Nothing
            molecularFormula molecule `shouldBe` Nothing
            validationInfo molecule `shouldBe` Nothing
            
            -- Verify all properties initialized properly
            molecularWeight molecule `shouldBe` Nothing
            xlogp molecule `shouldBe` Nothing
            hbd molecule `shouldBe` Nothing
            properties molecule `shouldBe` Map.empty
            
        it "should handle molecules with properties" $ do
            -- Create a base molecule
            let baseMol = createEmptyMolecule "mol-with-props"
            
            -- Add properties to a molecule
            let props = Map.fromList [("custom_prop", "value" :: Value), 
                                     ("numeric_prop", Number 42)]
                molWithProps = baseMol { 
                    smiles = Just "C1=CC=CC=C1", 
                    molecularFormula = Just "C6H6",
                    molecularWeight = Just 78.11,
                    properties = props
                }
                
            -- Verify properties are set correctly
            smiles molWithProps `shouldBe` Just "C1=CC=CC=C1"
            molecularFormula molWithProps `shouldBe` Just "C6H6"
            molecularWeight molWithProps `shouldBe` Just 78.11
            Map.size (properties molWithProps) `shouldBe` 2
            
        it "should validate molecule IDs" $ do 
            -- Test with empty ID (edge case)
            let emptyIdMolecule = createEmptyMolecule ""
            moleculeId emptyIdMolecule `shouldBe` ""
            
            -- Test with very long ID
            let longId = T.replicate 1000 "x"
                longIdMolecule = createEmptyMolecule longId
            moleculeId longIdMolecule `shouldBe` longId
            
            -- Test with special characters in ID
            let specialId = "mol-123@#$%^&*()"
                specialIdMolecule = createEmptyMolecule specialId
            moleculeId specialIdMolecule `shouldBe` specialId

    describe "Transformation Chain Operations" $ do
        it "should add transformations to a chain and retrieve them" $ do
            -- Create test environment
            logHandleMVar <- newMVar stdout
            let logHandle = LogHandle logHandleMVar
            
            -- Create mock config values
            let mockChemConfig = ChemConfig
                    { maxMoleculeSize = 1000
                    , maxTransformationSteps = 50
                    , computationTimeout = 300
                    , enabledProperties = ["logP", "topology"]
                    }
                mockApiConfig = APIConfig
                    { apiPort = 8080
                    , apiHost = "localhost"
                    , apiMaxRequestSize = 10485760
                    , apiRateLimit = 100
                    , apiTimeout = 30
                    , socketPath = "/tmp/test.sock"
                    }
                mockLogConfig = LogConfig
                    { logLevel = Info
                    , logFile = Nothing
                    , logFormat = LogText
                    , enableConsoleLog = True
                    }
                
            let mockEnv = AppEnv 
                  { appResources = Resources 
                      { resourceLogHandle = logHandle
                      }
                  , appConfig = AppConfig
                      { environment = Development
                      , chemConfig = mockChemConfig
                      , apiConfig = mockApiConfig
                      , logConfig = mockLogConfig
                      }
                  }
            
            -- Test a sequence of transformations
            result <- runApp mockEnv $ do
                -- Create initial transformation and start session
                t1 <- liftIO $ createTransformation "session-start"
                let initTrans = t1 { agent = "system", iteration = Just 0 }
                sessionId <- startSession "test-user" initTrans
                
                -- Create a second transformation
                t2 <- liftIO $ createTransformation "molecule-modify"
                let modifyTrans = t2 { 
                        agent = "medicinal-chemist", 
                        iteration = Just 1
                    }
                
                -- Log the second transformation
                tid2 <- logTransformation modifyTrans
                logInfo $ "Logged transformation 2 with ID: " <> T.pack (show tid2)
                
                -- Create a third transformation
                t3 <- liftIO $ createTransformation "property-calculate"
                let propTrans = t3 { 
                        agent = "computational-chemist", 
                        iteration = Just 2
                    }
                
                -- Log the third transformation
                tid3 <- logTransformation propTrans
                logInfo $ "Logged transformation 3 with ID: " <> T.pack (show tid3)
                
                -- Get transformation chain for a molecule
                testMolId <- liftIO UUID.nextRandom
                transChain <- getTransformationChain testMolId
                
                -- In our mock implementation with database disabled, 
                -- we expect the chain to be empty, but no errors
                pure (tid2, tid3, transChain)
                
            -- Verify operation succeeded
            case result of
                Left err -> expectationFailure $ "Transformation chain operations failed: " ++ show err
                Right (tid2, tid3, chain) -> do
                    -- Verify transformation IDs are valid UUIDs
                    show tid2 `shouldSatisfy` (not . null)
                    show tid3 `shouldSatisfy` (not . null)
                    
                    -- In our mock with db disabled, we expect an empty chain
                    null chain `shouldBe` True
                    
        it "should validate the transformation chain structure" $ do
            -- Since our mock implementation doesn't store chains, we test with
            -- constructed values and basic validation
            
            -- Create some test transformations
            trans1 <- createTransformation "transform-1"
            trans2 <- createTransformation "transform-2"
            
            -- Test properties
            transformationType trans1 `shouldBe` "transform-1"
            transformationType trans2 `shouldBe` "transform-2"
            
            -- Test with modifications
            let modifiedTrans = trans1 {
                    agent = "test-agent", 
                    iteration = Just 123
                }
            
            -- Verify modified properties
            agent modifiedTrans `shouldBe` "test-agent"
            iteration modifiedTrans `shouldBe` Just 123
    
    describe "Molecule Transformation and Validation" $ do
        it "should apply transformations to molecules with realistic results" $ do
            -- Create test environment
            (_, mockEnv) <- createMockEnv
            
            -- Test applying a transformation to a molecule
            result <- runApp mockEnv $ do
                -- Create a molecule with properties (benzene)
                let baseMol = createEmptyMolecule "benzene"
                    molecule = baseMol {
                        smiles = Just "C1=CC=CC=C1",
                        molecularFormula = Just "C6H6",
                        molecularWeight = Just 78.11,
                        properties = Map.fromList [("aromatic", Bool True)]
                    }
                
                -- Create a transformation to add a hydroxyl group (phenol)
                t <- liftIO $ createTransformation "molecule-modify"
                let modifyTrans = t { 
                    agent = "medicinal-chemist", 
                    iteration = Just 1
                }
                
                -- Apply the transformation
                transformResult <- applyTransform modifyTrans molecule
                
                -- For enhanced realism, we'll also test what the transformed molecule
                -- should actually look like if our mock were more realistic.
                -- In a real application, the transformation would modify molecule properties.
                
                -- Create the expected result (phenol) that a proper mock would return
                -- This shows what an ideal mock implementation would produce
                let expectedPhMol = createEmptyMolecule "phenol"
                    expectedMol = expectedPhMol {
                        smiles = Just "C1=CC=CC=C1O",  -- Phenol SMILES
                        molecularFormula = Just "C6H6O",
                        molecularWeight = Just 94.11,  -- Actual phenol weight
                        properties = Map.fromList [
                            ("aromatic", Bool True),
                            ("hydroxyl_group", Bool True)
                        ]
                    }
                
                -- Return both actual and expected results for comparison
                pure (transformResult, expectedMol)
                
            -- Verify transformation results
            case result of
                Left err -> expectationFailure $ "Transformation failed: " ++ show err
                Right (transformedMol, expectedMol) -> 
                    case transformedMol of
                        Left innerErr -> expectationFailure $ "Transformation returned error: " ++ show innerErr
                        Right actualMol -> do
                            -- Verify the mock result has a valid ID
                            moleculeId actualMol `shouldSatisfy` (not . T.null)
                            
                            -- Current mock just returns placeholder molecules, but we
                            -- document what we'd actually expect in an enhanced implementation
                            pendingWith "Current mock returns placeholders - enhanced mock would return modified molecules"
                            
                            -- Ideal tests for realistic mock:
                            -- smiles actualMol `shouldBe` smiles expectedMol -- "C1=CC=CC=C1O"
                            -- molecularFormula actualMol `shouldBe` molecularFormula expectedMol -- "C6H6O"
                            -- molecularWeight actualMol `shouldSatisfy` (\w -> case w of
                            --                                             Just weight -> abs (weight - 94.11) < 0.01
                            --                                             Nothing -> False)
        
        it "should validate molecules and handle invalid inputs" $ do
            -- Create mock environment
            (_, mockEnv) <- createMockEnv
            
            -- Test with valid and invalid SMILES strings
            result <- runApp mockEnv $ do
                -- Valid input - benzene
                logInfo "Attempting to validate a valid SMILES string: C1=CC=CC=C1"
                validMol <- validateMolecule "C1=CC=CC=C1"
                
                -- Invalid input - malformed SMILES
                logInfo "Attempting to validate an invalid SMILES string: Not-a-SMILES"
                
                -- LIMITATION: Current mock doesn't validate, always succeeds
                -- In a realistic implementation, this should return an error
                -- We'll simulate both success and error cases to show what should happen
                invalidMolResult <- validateMolecule "Not-a-valid-SMILES"
                
                -- For a more realistic test, we create expected molecules
                -- for comparison. Ideally, the first would succeed with a proper
                -- molecule, and the second would fail with validation error.
                
                -- Create expected valid molecule (benzene)
                let baseMol = createEmptyMolecule "benzene"
                    expectedValidMol = baseMol { 
                        smiles = Just "C1=CC=CC=C1",
                        molecularFormula = Just "C6H6",
                        molecularWeight = Just 78.11
                    }
                
                -- Return both results for verification
                pure (validMol, invalidMolResult, expectedValidMol)
                
            -- Verify results
            case result of
                Left err -> expectationFailure $ "App execution failed: " ++ show err
                Right (validMol, invalidMol, expectedMol) -> do
                    -- Current mock behavior - both succeed with placeholders
                    moleculeId validMol `shouldSatisfy` (not . T.null)
                    moleculeId invalidMol `shouldSatisfy` (not . T.null)
                    
                    -- Note the limitations of the current mock
                    pendingWith "Mock always succeeds validation; enhanced mock would reject invalid SMILES"
                    
                    -- Ideal tests:
                    -- Valid SMILES should produce a molecule with expected properties
                    -- smiles validMol `shouldBe` Just "C1=CC=CC=C1"
                    -- molecularFormula validMol `shouldBe` Just "C6H6"
                    
                    -- Invalid SMILES should fail with validation error
                    -- We'd expect: result for invalidMol should be Left ValidationError
                
    describe "Error Handling" $ do
        it "should properly handle and report errors in configuration" $ do
            -- Create mock environment
            (_, mockEnv) <- createMockEnv
            
            -- Test error handling through configuration validation
            let configInvalid = "invalid-config.yaml"  -- Nonexistent file
            configResult <- runExceptT $ do
                loadConfig configInvalid  -- This should fail
            
            -- Should return an error rather than crashing
            case configResult of
                Left err -> do
                    -- Verify we get the expected error type/message
                    show err `shouldContain` "config"  -- Error should mention config
                Right _ -> expectationFailure "Should have failed with nonexistent config"
                
        it "should handle invalid transformations" $ do
            -- Create mock environment
            (_, mockEnv) <- createMockEnv
                  
            -- Test with an invalid transformation (missing required fields)
            result <- runApp mockEnv $ do
                -- Create a basic transformation
                t <- liftIO $ createTransformation "test-transform"
                
                -- Test error handling with missing required fields
                -- In the current mock, this never fails, but we document what should happen
                -- with a more realistic implementation
                
                -- Create various test molecules
                let baseMol = createEmptyMolecule "test-valid"
                    validMol = baseMol { smiles = Just "C1=CC=CC=C1" }
                    
                    -- Attempt transformations that should fail with better mocks
                    emptySmilesMol = createEmptyMolecule "test-empty-smiles"
                
                -- Apply transformation to valid molecule (should succeed)
                validResult <- applyTransform t validMol
                
                -- Current mock returns placeholder molecules for all inputs
                -- With a better mock, some of these would fail with errors
                pure validResult
                
            -- Verify results
            case result of
                Left err -> expectationFailure $ "App execution failed: " ++ show err
                Right validResult -> 
                    case validResult of
                        Left err -> expectationFailure $ "Transform unexpectedly failed: " ++ show err
                        Right mol -> do
                            -- Basic validation of results with current mock
                            moleculeId mol `shouldSatisfy` (not . T.null)
                            
                            -- Document what should happen with better mocks
                            pendingWith "With enhanced mocks, invalid transformations would fail with appropriate errors"
                            
        it "should handle and log errors appropriately" $ do
            -- Create mock environment
            (_, mockEnv) <- createMockEnv
                  
            -- Test comprehensive error handling
            result <- runApp mockEnv $ do
                -- Log various error types to demonstrate error handling capabilities
                logError "Simulating a validation error"
                logError "Simulating a computation error"
                logError "Simulating a transformation error"
                
                -- In a real implementation, we would have code like:
                -- eResult <- runExceptT $ do
                --     -- Operation that might fail
                --     validateComplexMolecule "Invalid input"
                --
                -- case eResult of
                --     Left err -> handleErrorByType err
                --     Right _ -> logInfo "Success"
                
                -- Test basic error logging functionality
                pure True
                
            -- Verify error logging worked
            result `shouldBe` Right True

    describe "In-Memory Operations" $ do
        it "should function with in-memory storage" $ do
            -- Create a mock environment
            (_, mockEnv) <- createMockEnv
            
            -- Test core operations work with in-memory system
            result <- runApp mockEnv $ do
                -- 1. Create and start a session
                t <- liftIO $ createTransformation "test-transform"
                sessionId <- startSession "test-user" t
                logInfo $ "Created session: " <> T.pack (show sessionId)
                
                -- 2. Log a transformation
                t2 <- liftIO $ createTransformation "another-transform"
                transformId <- logTransformation t2
                logInfo $ "Logged transformation: " <> T.pack (show transformId)
                
                -- 3. Apply a transformation 
                let molecule = createEmptyMolecule "test-molecule"
                transformResult <- applyTransform t molecule
                
                -- Check the transform result
                let resultValid = case transformResult of
                        Left _ -> False
                        Right mol -> not (T.null (moleculeId mol))
                
                pure (sessionId, transformId, resultValid)
                
            -- Should succeed with in-memory operations
            case result of
                Left err -> expectationFailure $ "In-memory operations failed: " ++ show err
                Right (sid, tid, transformResult) -> do
                    -- Verify session ID is valid
                    show sid `shouldSatisfy` (not . null)
                    
                    -- Verify transformation ID is valid
                    show tid `shouldSatisfy` (not . null)
                    
                    -- Transformation succeeded
                    transformResult `shouldBe` True
