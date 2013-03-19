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

import MPFA

getMPFA nseg globalSd tsamps =   
  let mnvars = mpfa nseg tsamps

      inimpfa = L.fromList $ [50, 0.04, 0.1, (foldl1 max $ map fst mnvars)*1.1]
      laprs@(mn1,_,_) = laplaceApprox defaultAM {nmTol = 0.001} 
                                      (likeMPFA globalSd mnvars) [] [] 
                                      $ inimpfa
  in laprs

ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)

getMPFAN sess = do
  (t0s'::[Double], amps'::[Double],sds::[Double]) <- fmap ffile  $ readFile (sess++"/epsps3")
  (wf, wfAmp, sigs) <- getWf sess

  let tsamps = filter (getFilter sess) $ zip t0s' $ map (*wfAmp) amps'
  let mnvars = mpfa 200 tsamps
  let inimpfa = L.fromList $ [60, log 0.04, log 0.1, 1]
  let globalSd = runStat  meanF sds

  --fail $ show $ likeMPFA globalSd mnvars inimpfa

  let laprs@(mn1,_,_) = laplaceApprox defaultAM {nmTol = 0.0001} 
                                      (likeMPFA globalSd mnvars) [] [] 
                                      $ inimpfa
  return $ mn1 @> 0


main = do
  sess <- getSess "f5sm10"

  h <- openFile ("Figure5.tex") WriteMode 
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

  puts $ "Figure 5: session "++sess

  (wf, wfAmp, sigs) <- getWf sess

  --print wfAmp

  
  (t0s'::[Double], amps'::[Double],sds::[Double]) <- fmap ffile  $ readFile (sess++"/epsps3")
  let tsamps = filter (getFilter sess) $ zip t0s' $ map (*wfAmp) amps'
      t0s = map fst tsamps


  let weighCurve' = map (weighRegression tsamps ) t0s
      maxPcurve = foldl1 max weighCurve'
      pcurve = map (/(maxPcurve)) weighCurve'
  let globalSd = runStat  meanF sds

  --plotIt "pcurve" $ zip t0s pcurve :+: tsamps

  --plotIt "wf" wf

  --plotIt "wfs" $ concat [take 5 sigs, take 5 $ reverse sigs]

  vsamples::[L.Vector Double] <- fmap (drop (getBurnIn "sess") . map L.fromList . catMaybes . map safeRead . lines) $ readFile (take 6 sess++"/npq_samples")

  --print $ head vsamples
  --print $ last vsamples


  let ns = map (roundD . (@>0)) vsamples
      ps = map (unlogit . (@>2)) vsamples
      qs = map ((*wfAmp) . exp . (@>3)) vsamples
  {-plotIt "nplot" $ zip [(0::Double)..] $ ns
  plotIt "pplot" $ zip [(0::Double)..] $ ps
  plotIt "qplot" $ zip [(0::Double)..] $ qs

  plotIt "nhist" $ Histo 25 $ ns
  plotIt "phist" $ Histo 50 $ ps
  plotIt "qhist" $ Histo 50 $ qs

  plotIt "npcor" $ zip ns ps
  plotIt "nqcor" $ zip ns qs
  plotIt "pqcor" $ zip qs ps -}

  let mnvars = mpfa 200 tsamps

  --plotIt "mpfa" $ mnvars



  let inimpfa = L.fromList $ [60, log 0.04, log 0.1, 1]

  --fail $ show $ likeMPFA globalSd mnvars inimpfa

  let laprs@(mn1,_,_) = laplaceApprox defaultAM {nmTol = 0.0001} 
                                      (likeMPFA globalSd mnvars) [] [] 
                                      $ inimpfa

  --puts $ "likeMPFA init="++show (likeMPFA globalSd mnvars inimpfa)
  --puts $ "likeMPFA final="++show (likeMPFA globalSd mnvars mn1)

  --puts $ "\\begin{verbatim}"++show mn1++"\\end{verbatim}"

  --fail $ show laprs

  {-puts $ "max p = "++ show (foldl1 max $ L.toList $ L.subVector 4 (L.dim mn1 - 4) mn1) -}

  pcts <- forM [0..59] $ \i -> do
       putStrLn $ show i
       vsams::[L.Vector Double] <- fmap (drop (getBurnIn "sess") . map L.fromList . catMaybes . map safeRead . lines) $ readFile ("f6sm"++show i++"/npq_samples")
       --print $ head vsams
       --print $ last vsams
       let ns = map ((@>0)) vsams
           phis = map (unlogit . (@>2)) vsams
           pct =percentile' 50 ns
       print (pct, percentile' 0.8 phis)
       nmpfa <- getMPFAN $ "f6sm"++show i
       print nmpfa
       return $ ((pct, nmpfa), percentile' 0.8 phis, (HistoStyle "histeps" 25 $ ns, HistoStyle "histeps" 25 $ phis)) 
       --return $ HistoStyle "histeps" 25 $ ns 

  {-plotIt "cookn" $ zip [(0::Double)..] $ sort $ map fst3 pcts
  plotIt "cookp" $ zip [(0::Double)..] $ sort $ map snd3 pcts
  plotIt "cookhn" $ ManySup $ map ( fst . trd3)  pcts
  plotIt "cookhp" $ ManySup $ map ( snd . trd3)  pcts -}
  --print $ percentile' 0.3 [0.0, 0.001..1]

  let lns = Lines [LineWidth 1, LineType 1, LineColor "black"]  

  plotIt "fig5a" $ ((AxisLabels "time (s)" "Vm (mV)" $ lns $ sigsAlignZero $ concat[take 5 sigs, take 5 $ reverse sigs]) :||:
                      (AxisLabels "time (s)" "EPSP amplitude (mV)" tsamps)) :==: 
                      ((AxisLabels "N" "P" (zip ns ps)) :||:
                        (AxisLabels "N" "Q (mV?)") (zip ns qs))

  plotIt "fig5b" $ (AxisLabels "Q" "P" (zip qs ps) :||: AxisLabels "Variance" "Mean" mnvars) :==: 
                      (( Lines [LineWidth 1, LineType 1] $ ManySup $ map ( fst . trd3)  pcts) :||: (AxisLabels "Simulation number" "Percentile" $ zip [(0::Double)..] $ sort $ map (fst . fst3) pcts))


  plotIt "mpfahist" $ Histo 20 $ filter (<80) $ map (snd . fst3) pcts
  --psSam <- sampleIO $ fakesamP 50 2000
  --plotIt "pscurve" $ zip [(0::Double) ..] (psSam::[Double])
  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure5.tex"
  return ()

unlogit x= 1/(1+exp(-x))

safeRead x = case readsPrec 5 x of 
               [] -> Nothing
               (x,_):_ -> Just x

minus1 x = x+1

percentile x sorted_xs = 
  let n = realToFrac $ length sorted_xs
      go y [] = 0
      go y (x:xs) | y >= x    = realToFrac $ length xs + 1 
                  | otherwise = go y xs
  in go x sorted_xs / n

percentile' :: Double -> [Double] -> Double
percentile' x unsorted_xs = go x unsorted_xs 0 0 where
      go _ []     above below = below / (above+below)
      go y (x:xs) above below | y >= x    = go y xs above (below+1)
                              | otherwise = go y xs (above+1) below


roundD :: Double -> Double
roundD = realToFrac . round

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