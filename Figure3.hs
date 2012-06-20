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




{-measTrad sess = do
  x<- doesFileExist (sess++"/tradAmps")
  if x 
     then fmap read $ readFile (sess++"/tradAmps")
     else measTrad1 sess -}

measTrad sess = do
 cacheIn (sess++"/tradAmps") $ fmap catMaybes $ inEverySession $ whenContinues sess $ do
     rebaseRelativeTo sess
     vm <- signalsDirect "vm"
     sessionIdentifier <- getSessionName
     sessionStart <- getSessionStart
     spike <- events "spike" ()
     running <- durations "running" ()
     exclude <- durations "exclude" ()
     let swings = (\(lo,hi) -> abs(hi-lo)) <$$> sigStat (minF `both` maxF) vm
     let noGood = contains ((>5)//swings) running
     let spikeg = sortBy ( comparing (fst)) $ minInterval 0.1 $ notDuring exclude $ notDuring noGood spike
     let noiseSigs = take 50 $ limitSigs' (-0.11) (-0.01) $ around (spikeg) $ vm
     let epspSigs = during (durAroundEvent (0.03) 0.07 spikeg) vm 
     let aroundSpike = baseline (-0.003) 0.003 $ limitSigs' (-0.05) 0.05 $ around (spikeg) $ vm
     let ampPeak = snd $ head $ peak $ take 1 $ QU.averageSigs $ take 100 $ aroundSpike
     let tpeak = fst $ head $ peak $ take 1 $ QU.averageSigs $ aroundSpike
     let measDur  = measureBl (-0.003, 0.003) (tpeak-0.001,0.001+tpeak) vm spikeg
     --lift $ print (sessionIdentifier, tpeak)
     lift $ print (sessionIdentifier, length spikeg, length measDur)
     return $ Just measDur

main = do
  --let sess = "57246a"
  let sess = "22b152"
  h <- openFile ("Figure3.tex") WriteMode 
  let puts s = hPutStrLn  h $ s ++ "\n"
      plotIt nm obj = do gnuplotToPS (nm++".eps") $ obj
                         system $ "epstopdf "++nm++".eps"
                         puts $"\\includegraphics[width=16cm]{"++nm++"}\n\n"
 
  puts $ unlines     ["\\documentclass[11pt]{article}",
     "%include lhs2TeX.fmt",
     "%include polycode.fmt",
     "\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm]{geometry}",
     "\\usepackage{graphicx}",
     "\\begin{document}",
     "\n"]

  meass <- measTrad sess 

  (wf, wfAmp, sigs) <- getWf sess

  --puts $ "wfamp= "++show wfAmp++"\n"


  let measPts = (map (\((t1,t2),v)-> (t1,v)) $ concat meass)::[(Double, Double)]

      tfilt = filter $ \(t,v) -> t<1000

  --print $ length $  measPts

  --plotIt ("epsps_"++ take 6 sess) $ measPts

  let ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)

  let ampsFromFile npars = do
        (t0s'::[Double], amps::[Double],sds::[Double]) <- fmap ffile  $ readFile (sess++"/epsps"++show npars)
        let tsamps = zip t0s' $ map (*wfAmp) amps
        return tsamps

  tsamps1 <- ampsFromFile 1
  tsamps2 <- ampsFromFile 2
  tsamps3 <- ampsFromFile 3

  let rngIt x = Points [PointType 7] $ XRange 10 4000 $ YRange (-0.3) 1.6 x

  plotIt "pcurve1" $ (YLabel "Estimated amplitude (mV)" $ ("Window", rngIt measPts)) :||: (YTics [] $ ("Least squares", rngIt tsamps1))
        :==: ((YLabel "Estimated amplitude (mV)" $ ("OU", rngIt tsamps2)) :||: (YTics [] $ ("OU + obs", rngIt tsamps3)))

 
--  let failNewT = filter (\(t,v) -> v<0.5 && t > 700 && t < 800) tsamps 
  --puts $ "failNew = "++show failNewT

--  let failOldT = filter (\(t,v) -> v<0.4 && t > 350 && t < 400) measPts 
  --puts $ "failOld = "++show failOldT


  --puts $ show $ sigInfo $ head sigs

--  let Just failOldS = find (\(Signal dt t0 _)-> t0+100*dt > (fst $ failOldT!!0) ) sigs
--  let Just failNewS = find (\(Signal dt t0 _)-> t0+100*dt > (fst $ failNewT!!0) ) sigs

--  plotIt "failold" $ ("trad. fail", [failOldS]) :+: ("new fail", [failNewS])


  
  --plotIt "wf" wf
  --plotIt "tst" GnuplotTest

--  let var1 = runStat meanF $ map (localVar measPts) measPts
--  let var2 = runStat meanF $ map (localVar tsamps) tsamps

--  puts $ show (var1, var2)

  var12 <- forM datasess $ \sess1 -> do
     vars <- varDiff sess1
     --puts $ take 6 sess1++ "\t"++intercalate "\t" (map show vars)
     return vars 

  let vars = map (runStat meanF) $ transpose $ filter allNonZero var12


  puts $ show vars

  plotIt "varplot" $ XTicLabel [("Window", 0), ("Least sq",1),("OU",2),("OU+obs",3),("OU+obs+flat",4)] $ YLabel "Variance" $ YRange 0 0.06 $ Boxes $ zip [(0::Double)..] vars
  --puts $ show $ map fst var12
  --puts $ show $ map snd var12

  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure3.tex"
  return ()

allNonZero = all (>1e-8) 

varDiff sess = do
  meass <- measTrad sess  

  (wf, wfAmp, sigs) <- getWf sess

  let measPts = (map (\((t1,t2),v)-> (t1,v)) $ concat meass)::[(Double, Double)]
  let varTrad = runStat meanF $ map (localVar measPts) measPts

  vars <- forM [1..4] $ \npars -> do
     let ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)
     ex <- doesFileExist (sess++"/epsps"++show npars)
     if ex
        then do (t0s'::[Double], amps::[Double],sds::[Double]) <- fmap ffile  $ readFile (sess++"/epsps"++show npars)
                let tsamps = zip t0s' $ map (*wfAmp) amps
                return $ runStat meanF $ map (localVar tsamps) tsamps
        else return 0
  return $ varTrad : vars



autoCorrSig :: Signal Double -> Signal Double
autoCorrSig (Signal dt t0 vec) = Signal dt t0 $ L.fromList $ autoCorr $ L.toList vec

--http://www.bearcave.com/misl/misl_tech/wavelets/stat/index.html
autoCorr :: [Double] -> [Double]
autoCorr xs = acfs where
  n = length xs
  mean = runStat meanF xs
  denom = sum $ map (\xi-> (xi - mean)^2) xs
  f xt xlag = (xt - mean) * (xlag - mean)
  acf lag = (sum $ zipWith f xs (drop lag xs)) / denom
  acfs = map acf [0..(n-1)]

avSigs :: [Signal Double] -> Signal Double
avSigs sigs@((Signal dt t0 _):_) = Signal dt t0 meanVec where
  meanVec = runStat meanF $ map (\(Signal _ _ v) -> v) sigs

sigInfo (Signal dt t0 v) = (dt,t0, L.dim v)