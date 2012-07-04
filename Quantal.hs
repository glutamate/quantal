{-# LANGUAGE DeriveDataTypeable, ScopedTypeVariables, NoMonomorphismRestriction, ViewPatterns, PackageImports #-}
module Main where

import System.Environment
import Database
import Query hiding (io) 
import QueryTypes
import QueryUtils hiding (averageSigs)
import qualified QueryUtils as QU
import Data.Maybe
import Data.List
import Control.Monad
import System.Directory
import "probably" Math.Probably.RandIO
import QuantalHelp
import "baysig" Baysig.Estimate.RTS
import Data.Binary
import qualified Numeric.LinearAlgebra as L
import "probably" Math.Probably.MCMC
import "probably" Math.Probably.Sampler
import "probably" Math.Probably.FoldingStats
import qualified Math.Probably.PDF as PDF
import "probably" Math.Probably.NelderMead

import System.IO
import Data.Ord
import System.Posix.Directory
import System.Cmd

import Graphics.Gnewplot.Exec
import Graphics.Gnewplot.Types
import Graphics.Gnewplot.Style
import Graphics.Gnewplot.Panels
import Graphics.Gnewplot.Instances
import Graphics.Gnewplot.Histogram

import Control.Monad.Trans
import qualified Data.Array.IArray as IA


import "probably" Math.Probably.IterLap

main = do 
  sessApprox:dowhat:rest <- getArgs
  let remove_fail "fail_resolve" = sessApprox
      remove_fail s = s
  sess <- fmap remove_fail $ resolveApproxSession sessionsRootDir sessApprox

  createDirectoryIfMissing False $ take 6 sess

  when ('1' `elem` dowhat) $ epspSigs sess
  when ('2' `elem` dowhat) $ measNoise 4 sess

  when ('6' `elem` dowhat) $ simulate sess rest
  when ('3' `elem` dowhat) $ measAmps 4 sess
  when ('4' `elem` dowhat) $ measNPQ sess
--  when ('5' `elem` dowhat) $ summary sess

--  when ('7' `elem` dowhat) $ measAmps1 sess

  when ('8' `elem` dowhat) $ simulateAll [25, 50, 300] [1000,2500]
  when ('9' `elem` dowhat) $ simulateAll [200, 100] [1000, 2500]
 
  when ('t' `elem` dowhat) $ testCovM sess
  when ('a' `elem` dowhat) $ measNoise 3 sess
  when ('b' `elem` dowhat) $ measNoise 2 sess
  when ('c' `elem` dowhat) $ measNoise 1 sess
  when ('d' `elem` dowhat) $ measAmps 1 sess
  when ('e' `elem` dowhat) $ measAmps 3 sess
  when ('f' `elem` dowhat) $ measAmps 2 sess
  when ('z' `elem` dowhat) $ measAll sess

  when ('y' `elem` dowhat) $ simulate4 sess rest

  return ()


measAll sess = do 
--  system $ "~/.cabal/bin/quantal "++sess++" 1"
  system $ "~/.cabal/bin/quantal "++sess++" a"
  system $ "~/.cabal/bin/quantal "++sess++" b"
  system $ "~/.cabal/bin/quantal "++sess++" c"
  system $ "~/.cabal/bin/quantal "++sess++" 2"
  system $ "~/.cabal/bin/quantal "++sess++" d"
  system $ "~/.cabal/bin/quantal "++sess++" e"
  system $ "~/.cabal/bin/quantal "++sess++" f"
  system $ "~/.cabal/bin/quantal "++sess++" 3"
  return ()

testCovM sess = runRIO $ do
  let npars = 4
  LoadSignals sigs' <- io $ decodeFile $ take 6 sess++"/sigs_"++take 6 sess++"_noise"
  let initialV = L.join $ map L.fromList [ take npars [ 2::Double, -2, -6, -9], [-60]] 
      sigs = take 10 sigs'
      var1 = runStat varF $ L.toList $ snd $ observe $ head sigs
  io$ print $ var1
  io $ print $ exp (-5)
  io$ print $ tmax/dt
  let sigpts = snd $ observe $ head sigs
  io $ print $ L.dim sigpts
  io $ print $ posteriorNoiseV npars sigs initialV

  {-let covM = mkCovM np (-2) (2::Double) (-6) (-9)
      row1 = L.toList $ head $ L.toRows $ fst covM
  print $ take 10 $ row1 
  --print $ drop (np-10) $ row1 

  print $ take 10 $ IA.elems $ snd covM -}
--  let (inv,lndet)=L.invlndet $ fst covM 
--  print $ snd $ L.invlndet $ fst covM 
  return () 

simulateAll nss ntrs = do 
  forM_ ntrs $ \ntr -> do
   forM_ nss $ \ns -> do
     forM_ [1..5] $ \run -> do
       let ntrs = pad $ reverse $ drop 3 $ reverse $ show ntr
           nsstr = take 2 $ show ns           
           sessnm = ntrs++nsstr++show run
           runIt = "./Quantal "++sessnm++" 6 "++show ns++" "++show ntr
       system runIt
       writeFile (sessnm++"/sessions") $ show [sessnm]
       system $ "./Quantal "++sessnm++" 34"


simulate4 sess (nstarst:nsims:_)  = 
  forM_ [0..(read nsims-1)] $ \i -> do
    simulate (sess++show (i+read nstarst)) ["100", "2000"]
  

simulate sess rest = runRIO $ do
  let n = read $ rest!!0
  let ntrials = read $ rest!!1
  io $ createDirectoryIfMissing False $ take 6 sess
  sigs <- fmap ( map stagger . zip [1..] ) $ sample $ fakesam n ntrials
  io $ encodeFile (take 6 sess++"/sigs_"++take 6 sess++"_epsps") $ LoadSignals sigs 
  io $ writeFile (take 6 sess++"/noisePars3") $ show hatPars
  io $ writeFile (take 6 sess++"/sessions") $ show [sess]
  io $ system $ "~/.cabal/bin/quantal "++sess++" e" 
  io $ system $ "~/.cabal/bin/quantal "++sess++" 4"
  return ()





epspSigs sess = do 
  nms <- fmap catMaybes $ inEverySession $ whenContinues sess $ do
     rebaseRelativeTo sess
     vm <- signalsDirect "vm"
     sessionIdentifier <- getSessionName
     sessionStart <- getSessionStart
     spike <- events "spike" ()
     running <- durations "running" ()
     exclude <- durations "exclude" ()
     let swings = (\(lo,hi) -> abs(hi-lo)) <$$> sigStat (minF `both` maxF) vm
     let noGood = contains ((>3)//swings) running 
     let spikeg = sortBy ( comparing (fst)) $ minInterval 0.1 $ notDuring exclude $ notDuring noGood spike
     let noiseSigs = take 50 $ limitSigs' (0 - tmax- 0.01) (-0.01) $ around (spikeg) $ vm
     let epspSigs = during (durAroundEvent (0.03) 0.07 spikeg) vm 
     --let slopes = sigStat (fmap fst regressF) epspSigs
     let aroundSpike = baseline (-0.003) 0.003 $ limitSigs' (-0.05) 0.05 $ around (spikeg) $ vm
     let ampPeak = snd $ head $ peak $ take 1 $ QU.averageSigs $ take 100 $ aroundSpike
     let tpeak = fst $ head $ peak $ take 1 $ QU.averageSigs $ aroundSpike
     let measDur  = measureBl (-0.003, 0.003) (tpeak-0.015,0.015+tpeak) vm spikeg
     ask $ SaveSignals (take 6 sess++"/sigs_"++ take 6 sessionIdentifier ++ "_noise") noiseSigs
     ask $ SaveSignals (take 6 sess++"/sigs_"++ take 6 sessionIdentifier ++ "_epsps") epspSigs
     liftIO $ gnuplotToPNG (take 6 sess++"/epsps_"++ take 6 sessionIdentifier ++ ".png") $ measDur
     return $ Just sessionIdentifier
  writeFile (take 6 sess++"/sessions") $ show nms

measNoise npars sess = runRIO $ do 
  LoadSignals sigs' <- io $ decodeFile $ take 6 sess++"/sigs_"++take 6 sess++"_noise"
  let sigs = take 10 sigs'
      initparv = if npars ==1 then [-5] else [2::Double, 2, -6, -10]
      initialV = L.join $ map L.fromList [ take npars initparv, [sigAv $ head sigs]]    
  
  io$ print $ tmax/dt
  let sigpts = snd $ observe $ head sigs
  io $ print $ L.dim sigpts
  io $ print initialV
  io $ print $ posteriorNoiseV npars sigs initialV
  --io $ print $ posteriorNoiseV sigs initial1
  --io $ print $ posteriorNoiseV sigs initial2
  --fail "foo" 
  let fixed = [] -- [((i,j),0) | i <- [npars..1+npars], j <- [npars..1+npars], i/=j]
  laout@(init2,mbcor,_) <- io $ {-cacheIn (take 6 sess++"/laplaceNoise"++show npars) -} return $ laplaceApprox defaultAM {nmTol = 1, verboseNM = True,
                                                        initw = (\n -> if n<=npars-1 then 0.02 else 0.02)} 
                                             (posteriorNoiseV npars sigs) [] fixed initialV
  io $ print laout
--  io $ writeFile (take 6 sess++"/noiseLaplace"++show npars) $ show (init2, mbcor)
  let pdf0 = posteriorNoiseV npars sigs init2
  io $ print $ pdf0
  iniampar <- if notKosher mbcor
                 then do io $ print "starting from fresh"
                         sample $ initialAdaMet 200 (\n -> if n<=npars-1 then 5e-4 else 1e-3)  (posteriorNoiseV npars sigs) init2
                 else  initAdaMetFromCov 400 (posteriorNoiseV npars sigs) init2 0
                                                                          (L.scale 0.1 $ fromJust mbcor) 
 --  iniampar <- sample $ initialAdaMet 100 (\n -> if n<=3 then 1e-3 else 1e-3)  (posteriorNoiseV sigs) init2 

  --io$ print $ iniampar
  --let iniampar = AMPar init2 init2 (fromJust mbcor) 2.4 (pdf0) 0 0

  --froampar <- runAndDiscard 1500 (show . ampPar) iniampar $ adaMet False (posteriorNoiseV npars sigs)
 
  --io$ print $ froampar

  vsamples <- runAdaMetRIO 5000 True (iniampar {scaleFactor = 1.5}) (posteriorNoiseV npars sigs) 
  let vsmn = L.toList$ L.subVector 0 npars $ runStat meanF vsamples
{-      vsamTup = case vsmn of 
                   [sigma, logtheta, logobs, logsmooth] -> show (sigma, logtheta, logobs, logsmooth)
                   [sigma, logtheta, logobs] -> show (sigma, logtheta, logobs)
                   [sigma, logtheta] -> show (sigma, logtheta) -}
  io $ writeFile (take 6 sess++"/noisePars"++show npars) $ show vsmn 
  io $ writeFile (take 6 sess++"/noise_samples"++show npars) $ show vsamples
  return ()



notKosher (Nothing) = True
notKosher (Just mat) = any (nanOrInf) $ L.toList $ L.flatten mat

measAmps npars sess = runRIO $ do
  noisepars <- fmap read $ io $ readFile (take 6 sess++"/noisePars"++show npars)
  let covM = mkCovM np noisepars 
  let invDetails = invlndetC covM
  {-nms <- fmap read $ io $ readFile (take 6 sess++"/sessions")
  sigs <- fmap concat $ forM nms $ \sessNm-> do 
            LoadSignals sigs <- io $ decodeFile $ take 6 sess++"/sigs_"++take 6 sessNm++"_epsps" 
            return sigs
  let wf = baselineSig 0.003 $ averageSigs $ sigs -}
  (wf, _, sigs) <- io $ getWf sess
  h<- io $ openFile (take 6 sess++"/epsps"++show npars) WriteMode 
  forM_ sigs $ \sig@(Signal dt t0 _) -> do
      let initialV = L.fromList [-60,1]
      case laplaceApprox defaultAM {nmTol = 0.001} (posteriorSigV wf invDetails sig) [] [] initialV of

         (v, Just cor, smplx) -> do
                let amp = v @> 1
                    sd = sqrt $ (L.@@>) cor (1,1)
                io $ putStrLn $ "by Laplace: "++ show (t0,amp,sd)

                io $ hPutStrLn h $ show (t0, amp,sd)

         _             -> do
                vsamples <- nmAdaMet defaultAM (posteriorSigV wf invDetails sig) [] [] initialV
                let (amp,sd) = both meanF stdDevF `runStat` map (@>1) vsamples
                io $ print (t0,amp,sd)
                io $ hPutStrLn h $ show (t0, amp,sd)
      --plotPts $ zip [0..] $ map (@>1) vsamples
      return ()
  io $ hClose h
  return ()

{-measAmps1 sess = runRIO $ do
  (logtheta, sigma, logobs) <- fmap read $ io $ readFile (take 6 sess++"/noisePars")
  let covM = fillM (np,np) $
              \(i,j)-> ((covOU (exp logtheta) (sigma::Double)) (toD i)) (toD j)+ifObs i j (exp logobs)
  let invDetails = invlndet covM

  (wf, wfAmp, sigs) <- io $ getWf sess

  let sig@(Signal dt t0 _) = head sigs
  let initialV = L.fromList [-60,1]
{-  iniampar <- sample $ initialAdaMet 100 5e-3 (posteriorSigV wf invDetails sig) initialV
  froampar <- runAndDiscard 2000 (show . ampPar) iniampar $ adaMet False (posteriorSigV wf invDetails sig)
  vsamples<- runAdaMetRIO 3000 True froampar (posteriorSigV wf invDetails sig)
  let (amp,sd) = both meanF stdDevF `runStat` map (@>1) vsamples-}
--  io $ print (amp,sd)
  nmasams <- nmAdaMet defaultAM (posteriorSigV wf invDetails sig) [] [] initialV
  let (ampnm,sdnm) = both meanF stdDevF `runStat` map (@>1) nmasams
--  io $ print (amp,sd)
  io $ print (ampnm,sdnm) -}


  --plotPts $ zip [0..] $ map (@>1) vsamples
  return ()
 
measNPQ sess = runRIO $ do
  let ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)
  (t0s'::[Double], amps'::[Double],sds) <- io $ fmap ffile  $ readFile (take 6 sess++"/epsps3")
  let tsamps =  filter (getFilter sess) $ zip t0s' amps'
      --tsamps = filter ((<3) . (\(t, amp)-> zscore tsamps' (t,amp))) tsamps'
      t0s = map fst tsamps
      amps = map snd tsamps
  let weighCurve' = map (weighRegression tsamps ) t0s
      maxPcurve = foldl1 max weighCurve'
      pcurve = map (/(maxPcurve)) weighCurve'
  let globalSd = runStat ( meanF) sds

--  let initialV = (L.fromList [100,-10,0.6,-3]:: L.Vector Double)
 
--  io $ putStr "init ="
--  io $ print $ initialV

  io $ putStr "globalsd ="
  io $ print $ globalSd 

  io $ putStr "C 50 1 ="
  io $ print $ (choose 50 25 ,binomialProb 50 0.1 40, exp $ binomialLogProb 50 0.1 40)

  let fastPDF n v = posteriorNPQV amps pcurve globalSd $ L.join [ L.fromList [realToFrac n],  v]-- set n
      startN = getStartN sess

--  let npq@(maxV, maxN, _, smplx) = fastNPQ fastPDF  70 $  L.fromList [-10,0.5,-3.7]
  let npq@(maxV, maxN, _, smplx) = fastNPQup fastPDF  35 $  L.fromList [-10,0.5,-3.7]

  
  let maxFullV = L.join [ L.fromList [realToFrac maxN], maxV]

  io $ print $ maxFullV
 
  io $ print $ posteriorNPQV amps pcurve globalSd $ maxFullV

  --fail "boo!"

  let nsam = 100000
      nfrozen = 1000000


  iniampar <- sample $ initialAdaMet 5000 (const 5e-3) (posteriorNPQV amps pcurve globalSd) maxFullV
  io $ putStr "inipar ="
  io $ print $ iniampar 

  froampar <- runAndDiscard nsam (showNPQV') (iniampar {scaleFactor = 1.5}) $ adaMet False  (posteriorNPQV amps pcurve globalSd)
  io $ putStr "frozenpar ="
  io $ print $ froampar
  
  iniampar1 <- initAdaMetFromCov 20000 (posteriorNPQV amps pcurve globalSd) (ampMean froampar) 0
                                                                            (ampCov froampar) 

  vsamples <- runAdaMetRIO nsam True (iniampar {scaleFactor = 2.4}) (posteriorNPQV amps pcurve globalSd) 

  io $ writeFile (take 6 sess++"/npq_samples") $ show vsamples
  let (mean,sd) =  (both meanF stdDevF) `runStat` vsamples 
  io $ putStrLn $ intercalate "\t" $ map (minWidth 8) $ words "n cv phi q"
  io $ putStrLn $ showNPQV $ mean
  io $ putStrLn $ showNPQV sd  
  return ()

  {-let nsam = 40000
      nfrozen = 20000
  io $ print $ posteriorNPQV amps pcurve globalSd initialV
  iniampar <- sample $ initialAdaMet 500 5e-3 (posteriorNPQV amps pcurve globalSd) maxFullV
  io $ putStr "inipar ="
  io $ print $ ampPar iniampar 
  froampar <- runAndDiscard nsam (showNPQV') iniampar $ adaMet False (posteriorNPQV amps pcurve globalSd)
  io $ putStr "frozenpar ="
  io $ print $ ampPar froampar
  vsamples <- runAdaMetRIO nfrozen True froampar (posteriorNPQV amps pcurve globalSd) 
  io $ writeFile (take 6 sess++"/npq_samples") $ show vsamples
  let (mean,sd) =  (both meanF stdDevF) `runStat` vsamples 
  io $ putStrLn $ intercalate "\t" $ map (minWidth 8) $ words "n cv phi q"
  io $ putStrLn $ showNPQV $ mean
  io $ putStrLn $ showNPQV sd 
  return ()  -}

