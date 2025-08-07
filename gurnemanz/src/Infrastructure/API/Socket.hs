{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE CPP #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Infrastructure.API.Socket
    ( -- * Socket Server
      runSocketServer
    , handleClient
    -- * Message Types
    , Request(..)
    , Response(..)
    , TransformRequest(..)
    , SessionId
    ) where

import Control.Monad (forever, forM_, forM, when, void, mplus)
import Data.Maybe (catMaybes, isJust, isNothing, fromJust, mapMaybe, fromMaybe)
import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Error.Class (MonadError, throwError)
import Control.Concurrent (forkIO)
import Control.Exception (try, SomeException)
import qualified Data.List.NonEmpty as NE
import qualified Network.Socket as N
import qualified Network.Socket.ByteString as NB
import Data.ByteString (ByteString)
import qualified Data.ByteString.Lazy as LBS
import Data.Word (Word32)
import Data.Bits (shiftR, (.&.), shiftL, (.|.))
import Data.Aeson (ToJSON(..), FromJSON(..), Value(..), encode, decode, Value, toJSON, object, (.=), (.:), (.:?), (.!=), withObject, withText, parseJSON)
import System.Directory (removeFile)
import System.IO.Error (catchIOError)
import GHC.Generics (Generic)
import qualified Data.UUID as UUID
import qualified Data.UUID.V4 as UUID.V4
import Data.Text.Encoding (decodeUtf8)
import Data.List.NonEmpty (NonEmpty(..), (<|), toList)
import qualified Data.Map as Map
import Control.Concurrent.MVar (MVar, newMVar, modifyMVar, readMVar)
import qualified Data.Text as T
import Data.Text (Text)
import qualified Data.Vector as V
import Domain.Base.Types (Valid)
import Data.Aeson.Types (Parser)
import qualified Data.Aeson.Types as AT
import qualified Data.Aeson.KeyMap as KM
import qualified Data.ByteString as BS
import Data.Time.Clock (UTCTime, getCurrentTime)
import Data.Time.Format.ISO8601 (iso8601Show)
import Data.String (fromString)


-- Application specific imports
import Common.Types (LogLevel(..), LogContext(..), Resources(..))
import Core.Base (LogContext(..), LogLevel(..), MonadLogger(..), TransformationId, SessionId, UserId)
import Core.Error (Error(..))
import Infrastructure.Config.Types (APIConfig(..), logConfig)
import App.Types (App, runApp, AppEnv(..))  -- Import AppEnv record fields
import Domain.Molecule.Types (createEmptyMolecule, ValidatedMolecule, Valid)
import qualified Domain.Transform.Operations as Trans
import Domain.Base.Graph (MoleculeGraph(..))
import Domain.Transform.Types
    ( Transformation(..)
    , TransformationChain(..)
    , ChainStatus(..)
    , createTransformation
    , MethodDetails(..)
    )
import qualified Domain.Session.Types as Session
import Domain.Molecule.Types (ValidatedMolecule, createEmptyMolecule)
import Common.Types (LogLevel(..))

-- | Session store
type SessionStore = Map.Map SessionId TransformationChain

-- | Helper function to convert Maybe MethodDetails to Maybe Value
just_method_details :: Maybe Value -> Maybe Value
just_method_details (Just _) = Just (object [])
just_method_details Nothing = Nothing

-- | Helper function to lookup a Text key in an Object
-- We use a simple workaround since KM.fromText isn't directly available
lookup_with_text :: Text -> KM.KeyMap Value -> Maybe Value
lookup_with_text _ _ = Just (object [])  -- Simplified approach

-- | Transform request data
data TransformRequest = TransformRequest
    { trType :: Text
    , trUserMessage :: Maybe Text
    , trMethodDetails :: Maybe Value
    , trAgent :: Maybe Text         -- Agent field
    , trRationale :: Maybe Text     -- Rationale field
    , trProteinSequence :: Maybe Text -- Protein sequence field
    , trProteinPath :: Maybe Text   -- Protein structure file path
    , trIteration :: Maybe Int      -- Iteration number field
    , trUserId :: Maybe Text        -- User ID field
    } deriving (Show, Generic)

instance FromJSON TransformRequest where
    parseJSON = withObject "TransformRequest" $ \v ->
        TransformRequest
            <$> v .: "type"
            <*> v .:? "userMessage"
            <*> v .:? "methodDetails"
            <*> v .:? "agent"
            <*> v .:? "rationale"
            <*> v .:? "proteinSequence"
            <*> v .:? "proteinPath"
            <*> v .:? "iteration"
            <*> v .:? "userId"

instance ToJSON TransformRequest where
    toJSON tr = object
        [ "type" .= trType tr
        , "userMessage" .= trUserMessage tr
        , "methodDetails" .= trMethodDetails tr
        , "agent" .= trAgent tr
        , "rationale" .= trRationale tr
        , "proteinSequence" .= trProteinSequence tr
        , "proteinPath" .= trProteinPath tr
        , "iteration" .= trIteration tr
        , "userId" .= trUserId tr
        ]

-- | API Request types
data Request
    = CreateSession (Maybe TransformRequest) (Maybe [Value]) -- (transformation data, initial molecules)
    | ValidateMolecule SessionId Value
    | ApplyTransform SessionId UUID.UUID TransformRequest [Value] -- SessionId, parentTransformationId, transformation, molecules
    | RollbackTransformation SessionId Value
    | GetChain SessionId (Maybe UUID.UUID) (Maybe UUID.UUID) Bool -- SessionId, optional transformationId, optional moleculeId, includeMolecules
    deriving (Show, Generic)

-- | API Response types
data Response
    = Success Value
    | Failure Error
    deriving (Show, Generic)

-- Simplified ParseJSON instances for Request and Response (Slim version)
instance FromJSON Request where
    parseJSON = withObject "Request" $ \v -> do
        tag <- v .: "tag" :: Parser Text
        contents <- v .:? "contents"
        case tag of
            "CreateSession" -> case contents of
                Just obj -> case obj of
                    Object o -> do
                        transformData <- o .:? "transformation" :: Parser (Maybe TransformRequest)
                        moleculesVal <- o .:? "molecules" :: Parser (Maybe Value)
                        let molecules = case moleculesVal of
                                Just (Array arr) -> Just $ V.toList arr
                                _ -> Nothing
                        pure $ CreateSession transformData molecules
                    _ -> pure $ CreateSession Nothing Nothing
                _ -> pure $ CreateSession Nothing Nothing
                
            "GetChain" -> case contents of
                Just obj -> case obj of
                    Object o -> do
                        sidStr <- o .: "sessionId" :: Parser Text
                        tid <- o .:? "transformationId" :: Parser (Maybe Text)
                        mid <- o .:? "moleculeId" :: Parser (Maybe Text)
                        includeMols <- o .:? "includeMolecules" .!= False :: Parser Bool
                        
                        case UUID.fromText sidStr of
                            Nothing -> fail "Invalid UUID format for session ID"
                            Just sid -> do
                                let transId = case tid of
                                      Just t -> UUID.fromText t
                                      Nothing -> Nothing
                                let molId = case mid of
                                      Just m -> UUID.fromText m
                                      Nothing -> Nothing
                                pure $ GetChain sid transId molId includeMols
                    _ -> fail "Invalid GetChain request format"
                _ -> fail "Missing contents for GetChain"
                
            "ApplyTransform" -> case contents of
                Just obj -> case obj of
                    Array arr | V.length arr >= 3 -> do
                        let sessionIdVal = arr V.! 0
                        let parentTransformIdVal = arr V.! 1
                        let transformReqVal = arr V.! 2
                        let molecules = if V.length arr > 3 then V.toList $ V.drop 3 arr else []
                        
                        case (sessionIdVal, parentTransformIdVal) of
                            (String sid, String ptid) -> do
                                case (UUID.fromText sid, UUID.fromText ptid) of
                                    (Just sessionId, Just parentTransId) -> do
                                        transformReq <- parseJSON transformReqVal
                                        pure $ ApplyTransform sessionId parentTransId transformReq molecules
                                    _ -> fail "Invalid UUID format in ApplyTransform"
                            _ -> fail "SessionId or parentTransformationId is not a string"
                    _ -> fail "Invalid ApplyTransform request format"
                _ -> fail "Missing contents for ApplyTransform"
                
            "ValidateMolecule" -> case contents of
                Just (Array arr) | V.length arr == 2 -> do
                    let sessionIdVal = arr V.! 0
                    let moleculeVal = arr V.! 1
                    case sessionIdVal of
                        String sid -> do
                            case UUID.fromText sid of
                                Just sessionId -> pure $ ValidateMolecule sessionId moleculeVal
                                Nothing -> fail "Invalid UUID format for session ID"
                        _ -> fail "Session ID must be a string"
                _ -> fail "Invalid format for ValidateMolecule"
                
            "RollbackTransformation" -> case contents of
                Just (Array arr) | V.length arr == 2 -> do
                    let sessionIdVal = arr V.! 0
                    let transformationVal = arr V.! 1
                    case sessionIdVal of
                        String sid -> do
                            case UUID.fromText sid of
                                Just sessionId -> pure $ RollbackTransformation sessionId transformationVal
                                Nothing -> fail "Invalid UUID format for session ID"
                        _ -> fail "Session ID must be a string"
                _ -> fail "Invalid format for RollbackTransformation"
                
            _ -> fail $ "Unsupported request tag in slim version: " <> T.unpack tag

instance ToJSON Request where
    toJSON (CreateSession _ _) = object ["tag" .= ("CreateSession" :: Text)]
    toJSON _ = object ["tag" .= ("Unknown" :: Text)]  -- Simplified

instance ToJSON Response where
    toJSON (Success val) = object ["tag" .= ("Success" :: Text), "contents" .= val]
    toJSON (Failure err) = object ["tag" .= ("Failure" :: Text), "contents" .= err]

instance FromJSON Response where
    parseJSON = withObject "Response" $ \v -> do
        tag <- v .: "tag" :: Parser Text
        contents <- v .: "contents"
        case tag of
            "Success" -> pure $ Success contents
            "Failure" -> case AT.fromJSON contents of
                AT.Success err -> pure $ Failure err
                AT.Error msg -> fail $ "Invalid error format: " <> msg
            _ -> fail $ "Unknown response tag: " <> T.unpack tag

-- | Helper function to convert MoleculeGraph to ValidatedMolecule
toValidatedMolecule :: MoleculeGraph -> ValidatedMolecule Valid
toValidatedMolecule _ = createEmptyMolecule "placeholder-id"

-- Helper functions for working with time
-- removed formatISO8601 since we're using hardcoded dates instead

-- | Default logging context
defaultCtx :: LogContext
defaultCtx = LogContext
    { contextLevel = Info
    , contextComponent = "API"
    , contextCorrelationId = Nothing
    , contextTags = []
    }

-- | Run the socket server
runSocketServer :: AppEnv -> APIConfig -> IO ()
runSocketServer env APIConfig{..} = do
    -- Create a shared session store for all clients
    sessionStore <- newMVar Map.empty
    putStrLn $ "[DEBUG] Created new session store"
    
    catchIOError (removeFile socketPath) (const $ pure ())

    sock <- N.socket N.AF_UNIX N.Stream N.defaultProtocol
    N.bind sock (N.SockAddrUnix socketPath)
    N.listen sock 5

    putStrLn "[INFO] Socket server started"
    forever $ do
        (conn, _) <- N.accept sock
        -- Pass the shared session store to each client handler
        putStrLn $ "[DEBUG] Client connected, passing shared session store"
        _ <- forkIO $ handleClient env sessionStore conn
        pure ()

-- | Handle client connections
handleClient :: AppEnv -> MVar SessionStore -> N.Socket -> IO ()
handleClient env store conn =
    let loop :: IO ()
        loop = do
            -- Read 4-byte length prefix
            lengthBytes <- NB.recv conn 4
            if BS.null lengthBytes
                then do
                    putStrLn "[DEBUG] Client disconnected"
                    N.close conn
                    return ()
                else if BS.length lengthBytes < 4 then do
                    putStrLn $ "[DEBUG] Incomplete length prefix received: " ++ show lengthBytes
                    let err = Failure $ ValidationError "Incomplete length prefix"
                    sendResponseWithLength conn err
                    loop
                else do
                    -- Decode the message length from the 4 bytes
                    let messageLength = 
                            (fromIntegral (BS.index lengthBytes 0) `shiftL` 24) .|.
                            (fromIntegral (BS.index lengthBytes 1) `shiftL` 16) .|.
                            (fromIntegral (BS.index lengthBytes 2) `shiftL` 8) .|.
                            (fromIntegral (BS.index lengthBytes 3))
                    
                    putStrLn $ "[DEBUG] Expecting message of length: " ++ show messageLength ++ " bytes"
                    
                    -- Read the complete message in chunks
                    fullMsg <- readExactBytes conn messageLength
                    
                    if BS.length fullMsg /= messageLength
                        then do
                            putStrLn $ "[DEBUG] Failed to read complete message. Expected " 
                                    ++ show messageLength ++ " bytes, got " 
                                    ++ show (BS.length fullMsg)
                            let err = Failure $ ValidationError "Incomplete message received"
                            sendResponseWithLength conn err
                            loop
                        else do
                            putStrLn $ "[DEBUG] Received complete message of " ++ show messageLength ++ " bytes"
                            -- Attempt to decode the complete message
                            case decode $ LBS.fromStrict fullMsg of
                                Nothing -> do
                                    putStrLn $ "[DEBUG] Failed to decode: " ++ show fullMsg
                                    let err = Failure $ ValidationError "Invalid request format"
                                    sendResponseWithLength conn err
                                    loop
                                Just req -> do
                                    putStrLn $ "[DEBUG] Decoded request: " ++ show req
                                    
                                    -- Dump session store for debugging before processing
                                    sessions <- readMVar store
                                    putStrLn $ "[DEBUG] Current session store has " ++ show (Map.size sessions) ++ " sessions"
                                    forM_ (Map.toList sessions) $ \(sid, chain) -> do
                                        putStrLn $ "[DEBUG] Session " ++ show sid ++ " has " ++ 
                                                show (length (toList $ _steps chain)) ++ " steps"
                                    
                                    resp <- handleRequest env store req
                                    putStrLn $ "[DEBUG] Sending response: " ++ show resp
                                    
                                    -- Send the response with length prefix
                                    sendResponseWithLength conn resp
                                    loop
    in loop

createEmptyChain :: NonEmpty Transformation -> TransformationChain
createEmptyChain steps = TransformationChain
    { _steps = steps
    , _status = Active
    , _molecules = Map.empty
    }

-- | Helper function to read exactly n bytes from a socket
readExactBytes :: N.Socket -> Int -> IO ByteString
readExactBytes sock n = do
    go sock n []
  where
    go :: N.Socket -> Int -> [ByteString] -> IO ByteString
    go _ 0 chunks = pure $ BS.concat (reverse chunks)
    go s remaining chunks = do
        chunk <- NB.recv s (min 4096 remaining)
        if BS.null chunk
            then pure $ BS.concat (reverse chunks)  -- Connection closed prematurely
            else go s (remaining - BS.length chunk) (chunk : chunks)

-- | Helper function to send a response with length prefix
sendResponseWithLength :: (ToJSON a) => N.Socket -> a -> IO ()
sendResponseWithLength conn resp = do
    let responseBytes = LBS.toStrict $ encode resp
    let responseLength = fromIntegral (BS.length responseBytes) :: Word32
    putStrLn $ "[DEBUG] Response size: " ++ show responseLength ++ " bytes"
    
    -- Send 4-byte length prefix as big-endian bytes
    let lengthBytes = BS.pack [ 
            fromIntegral (responseLength `shiftR` 24) .&. 0xFF,
            fromIntegral (responseLength `shiftR` 16) .&. 0xFF,
            fromIntegral (responseLength `shiftR` 8) .&. 0xFF,
            fromIntegral responseLength .&. 0xFF
          ]
    NB.sendAll conn lengthBytes
    
    -- Send actual response
    NB.sendAll conn responseBytes

-- | Parse method details from JSON
parseMethodDetails :: Value -> Maybe MethodDetails
parseMethodDetails val = decode $ encode val

-- | Create a simplified transformation
createSimpleTransformation :: Text -> IO Transformation
createSimpleTransformation transformType = do
    tid <- UUID.V4.nextRandom
    pure $ Transformation
        { transformationId = tid
        , transformationType = transformType
        , userMessage = Nothing
        , agentResponse = Nothing
        , methodDetails = Nothing
        , inputMoleculeIds = []
        , outputMoleculeIds = []
        , agent = "system"
        , proteinSequence = Nothing
        , proteinPath = Nothing
        , iteration = Nothing
        }

-- | Add a transformation to a chain
addTransformation :: TransformationChain -> Transformation -> TransformationChain
addTransformation chain newTransform = 
    -- Create a deep copy of the chain with the new transformation at the head
    let newStep = newTransform <| _steps chain
    in chain { _steps = newStep }

-- Debug information about the chain
addTransformationDebug :: TransformationChain -> Transformation -> IO TransformationChain
addTransformationDebug chain newTransform = do
    -- Log existing chain information
    let existingSteps = toList $ _steps chain
    putStrLn $ "[DEBUG] Current chain has " ++ show (length existingSteps) ++ " steps"
    
    -- Create a deep copy of the chain with the new transformation at the head
    let newStep = newTransform <| _steps chain
    let newChain = chain { _steps = newStep }
    
    -- Log updated chain
    let updatedSteps = toList $ _steps newChain
    putStrLn $ "[DEBUG] Updated chain now has " ++ show (length updatedSteps) ++ " steps"
    putStrLn $ "[DEBUG] First step ID: " ++ show (transformationId (NE.head $ _steps newChain))
    
    return newChain

-- | Handle incoming requests (core transformation operations)
handleRequest :: AppEnv -> MVar SessionStore -> Request -> IO Response
handleRequest env store = \case
    CreateSession transformReq molecules -> do
        sid <- liftIO UUID.V4.nextRandom
        liftIO $ putStrLn $ "[DEBUG] CreateSession: Creating new session with ID: " ++ show sid
        
        -- Create initial transformation (either from supplied data or default)
        initialTrans <- case transformReq of
            Just tr -> do
                tid <- liftIO UUID.V4.nextRandom
                let transformType = trType tr
                let userMsg = trUserMessage tr
                let methodDetails = parseMethodDetails =<< trMethodDetails tr
                let agentName = maybe "user" id (trAgent tr)
                let rationale = trRationale tr
                let userId = trUserId tr
                
                liftIO $ putStrLn $ "[DEBUG] CreateSession: Using supplied transformation data: " ++ show transformType
                liftIO $ putStrLn $ "[DEBUG] CreateSession: User ID: " ++ show userId
                
                pure $ Transformation
                    { transformationId = tid
                    , transformationType = transformType
                    , userMessage = userMsg
                    , agentResponse = Nothing
                    , methodDetails = methodDetails
                    , inputMoleculeIds = []
                    , outputMoleculeIds = []
                    , agent = agentName
                    , proteinSequence = trProteinSequence =<< transformReq
                    , proteinPath = trProteinPath =<< transformReq
                    , iteration = trIteration =<< transformReq
                    }
            Nothing -> do
                initialT <- liftIO $ createSimpleTransformation "init"
                liftIO $ putStrLn $ "[DEBUG] CreateSession: Using default transformation"
                pure initialT
        
        liftIO $ putStrLn $ "[DEBUG] CreateSession: Initial transformation ID: " ++ show (transformationId initialTrans)
        
        -- Process molecules if provided
        moleculeIds <- case molecules of
            Just mols -> do
                liftIO $ putStrLn $ "[DEBUG] CreateSession: Processing " ++ show (length mols) ++ " initial molecules"
                
                -- Generate molecule IDs and create molecules
                forM mols $ \molData -> do
                    molId <- liftIO UUID.V4.nextRandom
                    liftIO $ putStrLn $ "[DEBUG] CreateSession: Generated molecule ID: " ++ show molId
                    pure molId
                
            Nothing -> do
                liftIO $ putStrLn $ "[DEBUG] CreateSession: No initial molecules provided"
                pure []
                
        -- If we have molecules, update the initialTrans to include them as outputs
        let finalTrans = initialTrans { outputMoleculeIds = moleculeIds }
        let initialChain = createEmptyChain (finalTrans :| [])
        
        -- Debug the initial chain structure
        liftIO $ putStrLn $ "[DEBUG] CreateSession: Initial chain created with " ++ 
                           show (length (toList $ _steps initialChain)) ++ " steps"
        liftIO $ putStrLn $ "[DEBUG] CreateSession: Initial transformation has " ++ 
                           show (length (outputMoleculeIds finalTrans)) ++ " output molecules"
        
        -- All data stored in-memory
        -- Add the session to the store atomically
        response <- modifyMVar store $ \sessions -> do
            -- Check for existing sessions
            liftIO $ putStrLn $ "[DEBUG] CreateSession: Store currently has " ++ show (Map.size sessions) ++ " sessions"
            
            -- Insert the new session
            let newSessions = Map.insert sid initialChain sessions
            liftIO $ putStrLn $ "[DEBUG] CreateSession: Store now has " ++ show (Map.size newSessions) ++ " sessions"
            
            -- Prepare enhanced transformation with generated ID
            let enhancedTransformRequest = case transformReq of
                    Just tr -> object $
                        [ "type" .= trType tr
                        , "userMessage" .= trUserMessage tr
                        , "methodDetails" .= trMethodDetails tr
                        , "agent" .= maybe "user" id (trAgent tr)
                        , "rationale" .= trRationale tr
                        , "transformationId" .= transformationId finalTrans  -- Add the generated transformationId
                        , "userId" .= trUserId tr  -- Include the user ID in the response
                        ]
                    Nothing -> object 
                        [ "type" .= ("init" :: Text)
                        , "transformationId" .= transformationId finalTrans
                        ]
            
            -- Prepare enhanced molecules with generated IDs
            let enhancedMolecules = case molecules of
                    Just mols -> zipWith (\molData molId -> object [ 
                                "moleculeId" .= molId,
                                "originalData" .= molData 
                            ]) mols moleculeIds
                    Nothing -> []
            
            -- Create a response with sessionId and the request data with embedded IDs
            let responseValue = object
                    [ "sessionId" .= sid  -- Keep sessionId at the top level since it's not part of the request
                    , "transformation" .= enhancedTransformRequest
                    , "molecules" .= enhancedMolecules
                    ]
            
            pure (newSessions, Success responseValue)
            
        return response
    
    ApplyTransform sid parentTransformId transformReq molecules -> do
        -- Parse transformation details from the request
        let transformType = trType transformReq
            userMsg = trUserMessage transformReq
            methodDetails = parseMethodDetails =<< trMethodDetails transformReq
            agentName = maybe "system" id (trAgent transformReq)
            reason = trRationale transformReq
            proteinSeq = trProteinSequence transformReq
            proteinPath = trProteinPath transformReq
            iteration = trIteration transformReq
        
        -- Create a new transformation with a unique ID
        tid <- liftIO UUID.V4.nextRandom
        
        -- Debug logs to see what we're receiving
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Creating transformation with type: " ++ show transformType
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Agent: " ++ show agentName
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Session ID: " ++ show sid
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Parent transformation ID: " ++ show parentTransformId
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: New transformation ID: " ++ show tid
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Number of molecules: " ++ show (length molecules)
        
        -- Get resources
        let resources = appResources env
            logHandle = resourceLogHandle resources
            config = logConfig $ appConfig env
            ctx = LogContext Info "ApplyTransform" Nothing []
        
        -- Extract parent molecule IDs from the molecules data
        let extractParentMoleculeId :: Value -> Maybe UUID.UUID
            extractParentMoleculeId molData = case molData of
                Object o -> case AT.parseEither (\obj -> obj .:? "parentMoleculeId") o of
                    Right (Just idStr) -> UUID.fromText idStr
                    _ -> Nothing
                _ -> Nothing
        
        let parentMoleculeIds = mapMaybe extractParentMoleculeId molecules
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Parent molecule IDs: " ++ show parentMoleculeIds
        
        -- Verify session store before update
        beforeStore <- readMVar store
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Store before update has " ++ show (Map.size beforeStore) ++ " sessions"
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Session exists in store: " ++ show (Map.member sid beforeStore)
        
        -- Generate IDs for the new molecules
        newMoleculeIds <- forM molecules $ \_ -> liftIO UUID.V4.nextRandom
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Generated molecule IDs: " ++ show newMoleculeIds
        
        -- All molecules stored in-memory with the transformation chain
        -- Update the session store with the new transformation
        result <- modifyMVar store $ \sessions -> do
            liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Inside modifyMVar, store has " ++ show (Map.size sessions) ++ " sessions"
            
            case Map.lookup sid sessions of
                Nothing -> do
                    liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Session not found: " ++ show sid
                    pure (sessions, Failure $ ValidationError "Session not found")
                Just chain -> do
                    -- Debug current chain state
                    liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Before Apply: Chain has " ++ 
                                      show (length $ toList $ _steps chain) ++ " steps"
                    
                    -- Verify parent transformation exists
                    let allTransforms = toList $ _steps chain
                    let parentExists = any (\t -> transformationId t == parentTransformId) allTransforms
                    
                    if not parentExists
                        then do
                            liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Parent transformation not found: " ++ show parentTransformId
                            pure (sessions, Failure $ ValidationError "Parent transformation not found")
                        else do
                            -- Create new transformation with input molecule IDs from parent molecules
                            let transform = Transformation
                                    { transformationId = tid
                                    , transformationType = transformType
                                    , userMessage = userMsg
                                    , agentResponse = Nothing
                                    , methodDetails = methodDetails
                                    , inputMoleculeIds = parentMoleculeIds
                                    , outputMoleculeIds = newMoleculeIds
                                    , agent = agentName
                                    , proteinSequence = trProteinSequence transformReq
                                    , proteinPath = trProteinPath transformReq
                                    , iteration = trIteration transformReq
                                    }
                            
                            -- Use addTransformation to add to the chain
                            let newChain = Trans.addTransformation chain transform
                            
                            -- Create molecules with data to add to the chain
                            let moleculePairs = zipWith (\molData molId -> (molId, molData)) molecules newMoleculeIds
                            
                            -- Add molecules to the chain
                            let newChainWithMolecules = Trans.addMolecules newChain moleculePairs
                            
                            -- Debug the updated chain
                            liftIO $ putStrLn $ "[DEBUG] ApplyTransform: After Apply: Chain now has " ++ 
                                              show (length $ toList $ _steps newChainWithMolecules) ++ " steps"
                            liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Added " ++ 
                                              show (length moleculePairs) ++ " molecules to the chain"
                            
                            -- Log every transformation in the chain for debugging
                            let allSteps = toList $ _steps newChainWithMolecules
                            liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Chain has " ++ 
                                              show (Map.size $ _molecules newChainWithMolecules) ++ " molecules stored"
                            liftIO $ forM_ (zip [1..] allSteps) $ \(i, t) -> do
                                putStrLn $ "[DEBUG] ApplyTransform:   Step " ++ show i ++ ": " ++ 
                                          "id=" ++ show (transformationId t) ++ 
                                          ", type=" ++ show (transformationType t) ++ 
                                          ", agent=" ++ show (agent t) ++
                                          ", in=" ++ show (inputMoleculeIds t) ++
                                          ", out=" ++ show (outputMoleculeIds t)
                            
                            -- Create an enhanced response with the original data and added IDs
                            -- Add original transform request but with newly generated IDs
                            let enhancedTransformRequest = object $
                                    [ "type" .= transformType
                                    , "userMessage" .= userMsg
                                    , "methodDetails" .= trMethodDetails transformReq
                                    , "agent" .= agentName
                                    , "rationale" .= reason
                                    , "transformationId" .= tid  -- Add the generated transformationId
                                    , "proteinSequence" .= trProteinSequence transformReq  -- Add the protein sequence
                                    , "proteinPath" .= trProteinPath transformReq  -- Add the protein structure file path
                                    , "iteration" .= trIteration transformReq  -- Add the iteration number
                                    ]
                            
                            -- Add molecule IDs to each molecule in the original data
                            let enhancedMolecules = zipWith (\molData molId -> object [ 
                                        "moleculeId" .= molId,
                                        "originalData" .= molData 
                                    ]) molecules newMoleculeIds
                            
                            -- Build the response with just the original request data with IDs embedded
                            let responseObject = object
                                    [ "transformation" .= enhancedTransformRequest
                                    , "molecules" .= enhancedMolecules
                                    ]
                            
                            -- Insert updated chain into session store
                            let updatedSessions = Map.insert sid newChainWithMolecules sessions
                            
                            -- Verify the chain was actually updated in the map
                            case Map.lookup sid updatedSessions of
                                Nothing -> liftIO $ putStrLn $ "[DEBUG] ApplyTransform: ERROR - Session missing after update"
                                Just updatedChain ->
                                    liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Session store updated, chain steps: " ++ 
                                                     show (length $ toList $ _steps updatedChain) ++ ", molecules: " ++
                                                     show (Map.size $ _molecules updatedChain)
                            
                            pure (updatedSessions, Success responseObject)
        
        -- Verify session store after update
        afterStore <- readMVar store
        liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Store after update has " ++ show (Map.size afterStore) ++ " sessions"
        
        case Map.lookup sid afterStore of
            Nothing -> liftIO $ putStrLn $ "[DEBUG] ApplyTransform: ERROR - Session not found after update"
            Just finalChain -> do
                liftIO $ putStrLn $ "[DEBUG] ApplyTransform: Final chain has " ++ 
                                  show (length $ toList $ _steps finalChain) ++ " steps"
                liftIO $ putStrLn $ "[DEBUG] ApplyTransform: First step ID: " ++ 
                                  show (transformationId $ NE.head $ _steps finalChain)
        
        pure result
    
    GetChain sid mTransformId mMoleculeId includeMolecules -> do
        -- Retrieve the chain from the session store
        sessions <- readMVar store
        liftIO $ putStrLn $ "[DEBUG] GetChain: Looking for session " ++ show sid
        
        -- Debug all sessions in the store
        liftIO $ putStrLn $ "[DEBUG] GetChain: Current session store has " ++ show (Map.size sessions) ++ " sessions"
        liftIO $ putStrLn $ "[DEBUG] GetChain: Optional transformationId: " ++ show mTransformId
        liftIO $ putStrLn $ "[DEBUG] GetChain: Optional moleculeId: " ++ show mMoleculeId
        liftIO $ putStrLn $ "[DEBUG] GetChain: Include molecules: " ++ show includeMolecules
        
        -- Check if we have neither a transformationId nor a moleculeId
        when (isNothing mTransformId && isNothing mMoleculeId) $
            liftIO $ putStrLn $ "[DEBUG] GetChain: No transformationId or moleculeId specified, returning full chain"
        
        -- List all sessions for debugging
        liftIO $ putStrLn $ "[DEBUG] GetChain: All session IDs in store: " ++ show (Map.keys sessions)
        
        -- Check if session exists
        case Map.lookup sid sessions of
            Nothing -> do
                liftIO $ putStrLn $ "[DEBUG] GetChain: Session " ++ show sid ++ " not found"
                pure $ Failure $ ValidationError "Session not found"
            Just chain -> do
                -- Get the raw steps from the chain
                let rawSteps = _steps chain
                let allSteps = toList rawSteps
                liftIO $ putStrLn $ "[DEBUG] GetChain: Raw chain has " ++ show (length allSteps) ++ " steps"
                
                -- Handle based on which path we're taking
                if isJust mMoleculeId
                    then do
                        -- Molecule Chain Path
                        let molId = fromJust mMoleculeId
                        liftIO $ putStrLn $ "[DEBUG] GetChain: Building molecule chain for: " ++ show molId
                        
                        -- Find transformations that produced this molecule
                        let producingTransforms = filter (\t -> molId `elem` outputMoleculeIds t) allSteps
                        
                        if null producingTransforms
                            then do
                                liftIO $ putStrLn $ "[DEBUG] GetChain: No transformations found for molecule: " ++ show molId
                                pure $ Failure $ ValidationError "No transformations found for molecule"
                            else do
                                let transform = head producingTransforms
                                liftIO $ putStrLn $ "[DEBUG] GetChain: Found transformation producing molecule: " ++ 
                                                  show (transformationId transform)
                                
                                -- Build the molecule chain by following parentMoleculeIds
                                -- Try to get the actual molecule data from our stored molecules
                                
                                -- Create molecule data from our stored data if available
                                let moleculeData = case Map.lookup molId (_molecules chain) of
                                        Just molData -> object
                                            [ "moleculeId" .= molId
                                            , "transformationId" .= transformationId transform
                                            , "data" .= molData  -- Include the stored molecule data
                                            ]
                                        Nothing -> object
                                            [ "moleculeId" .= molId
                                            , "transformationId" .= transformationId transform
                                            ]
                                
                                -- Construct response with chain and optionally include full molecule data
                                let chainResponse = if includeMolecules
                                        then object
                                            [ "chain" .= [moleculeData]
                                            , "molecules" .= Map.filterWithKey (\k _ -> k == molId) (_molecules chain)
                                            ]
                                        else object
                                            [ "chain" .= [moleculeData]
                                            ]
                                
                                liftIO $ putStrLn $ "[DEBUG] GetChain: Returning molecule chain with " ++
                                                  (if includeMolecules then "full molecule data" else "references only")
                                
                                pure $ Success chainResponse
                    
                    else if isJust mTransformId
                        then do
                            -- Transformation Chain Path - from root to specified transform
                            let transformId = fromJust mTransformId
                            liftIO $ putStrLn $ "[DEBUG] GetChain: Building transformation chain for: " ++ show transformId
                            
                            -- Find the specific transformation
                            let matchingTransform = filter (\t -> transformationId t == transformId) allSteps
                            
                            if null matchingTransform
                                then do
                                    liftIO $ putStrLn $ "[DEBUG] GetChain: Transformation not found: " ++ show transformId
                                    pure $ Failure $ ValidationError "Transformation not found"
                                else do
                                    -- In a real implementation, we would trace backward from here
                                    -- For now, we'll just return all transformations up to and including this one
                                    
                                    let transformIndex = length allSteps - 
                                                        length (dropWhile (\t -> transformationId t /= transformId) (reverse allSteps))
                                    let relevantTransforms = take (transformIndex + 1) (reverse allSteps)
                                    
                                    liftIO $ putStrLn $ "[DEBUG] GetChain: Found " ++ 
                                                      show (length relevantTransforms) ++ " transformations in chain"
                                    
                                    -- Build response with the chain and optionally include full molecule data
                                    let chainResponse = if includeMolecules
                                            then object
                                                [ "chain" .= relevantTransforms
                                                , "molecules" .= _molecules chain  -- Include full molecule data
                                                ]
                                            else object
                                                [ "chain" .= relevantTransforms
                                                ]
                                    
                                    pure $ Success chainResponse
                        
                        else do
                            -- Neither specified, return full chain
                            let transformations = reverse allSteps  -- Chronological order
                            
                            -- Debug each transformation in detail
                            forM_ (zip [1..] transformations) $ \(i, t) -> do
                                liftIO $ putStrLn $ "[DEBUG] GetChain:   Step " ++ show i ++ ": " ++ 
                                                  "id=" ++ show (transformationId t) ++ 
                                                  ", type=" ++ show (transformationType t) ++ 
                                                  ", agent=" ++ show (agent t) ++ 
                                                  ", in=" ++ show (inputMoleculeIds t) ++ 
                                                  ", out=" ++ show (outputMoleculeIds t)
                            
                            -- Build response with all steps and optionally include molecule data
                            let chainResponse = if includeMolecules
                                    then object
                                        [ "chain" .= transformations
                                        , "status" .= _status chain
                                        , "molecules" .= _molecules chain  -- Include full molecule data
                                        ]
                                    else object
                                        [ "chain" .= transformations
                                        , "status" .= _status chain
                                        ]
                            
                            liftIO $ putStrLn $ "[DEBUG] GetChain: Returning full chain with " ++ 
                                               show (length transformations) ++ " steps" ++
                                               (if includeMolecules then " and " ++ show (Map.size $ _molecules chain) ++ " molecules" else "")
                            
                            pure $ Success chainResponse
    
    RollbackTransformation sid transformId -> do
        -- Simple implementation to demonstrate rollback
        sessions <- readMVar store
        case Map.lookup sid sessions of
            Nothing -> pure $ Failure $ ValidationError "Session not found"
            Just chain -> 
                -- In a real implementation, we would filter the chain to keep only
                -- transforms up to and including the specified transformId
                pure $ Success $ toJSON chain

-- | Update a transformation in a chain
updateTransformationInChain :: TransformationChain -> UUID.UUID -> AT.Object -> TransformationChain
updateTransformationInChain chain tid updateData =
    let agentResponse = case AT.parseEither (.: "agentResponse") updateData of
            Right resp -> Just resp
            Left _ -> Nothing
            
        outputIds = case AT.parseEither (.: "outputMoleculeIds") updateData of
            Right ids -> ids
            Left _ -> []
            
        agent = case AT.parseEither (.: "agent") updateData of
            Right a -> a
            Left _ -> "system"
            
        -- Debug the transformation being updated    
        updateStep t = 
            if transformationId t == tid
            then 
                let updated = t { agentResponse = agentResponse 
                               , outputMoleculeIds = outputIds
                               , agent = agent
                               }
                in updated
            else t
            
        -- Apply the update to all steps in the chain
        steps' = fmap updateStep (_steps chain)
    in chain { _steps = steps' }