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

main = do
  let sess = "57246a"
  h <- openFile ("Figure1.tex") WriteMode 
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

  meass <- fmap catMaybes $ inEverySession $ whenContinues sess $ do
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
     --lift $ print (sessionIdentifier, length spikeg, length measDur)
     return $ Just measDur

  nms <- fmap read $  readFile (take 6 sess++"/sessions")
  sigs <- fmap concat $ forM nms $ \sessNm-> do 
            LoadSignals sigs <- decodeFile $ take 6 sess++"/sigs_"++take 6 sessNm++"_epsps" 
            return sigs
  let wf@(Signal _ _ sv) =  averageSigs $ sigs
  let wfAmp = foldl1' max $ L.toList sv
  puts $ "Figure 1\n\n"

  let plotA = AxisLabels "time (s)" "membrane voltage (mV)" $ 
                  (Lines [LineWidth 1, LineType 1, LineColor "red"] $ sigsAlignZero $ take 10 $ drop 250 sigs) :+: 
                  (Lines [LineWidth 4, LineType 1, LineColor "black"] $ sigAlignZero $ smap (\v -> v-2.5) wf)


--  plotIt "wf" wf


  let measPts = (map (\((t1,t2),v)-> (t1,v)) $ concat meass)::[(Double, Double)]

      tfilt = filter $ \(t,v) -> t<1000


  let plotB = AxisLabels "time (s)" "EPSP (mV)" $ XRange 0 10000 measPts

  let plotC = AxisLabels "mean (mV)" "variance (mV^2)" $ mpfa 100 measPts

  let plotD = AxisLabels "mean (mV)" "variance (mV^2)" $mpfa 500 measPts

  plotIt "fig1" $ (A plotA :||: B plotB) :==: (C plotC :||: D plotD)
  plotIt "gnuplot_test" $ GnuplotTest

--  plotIt ("epsps_"++ take 6 sess) $ measPts

  

  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure1.tex"
  return ()

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