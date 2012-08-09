{-# LANGUAGE DeriveDataTypeable, ScopedTypeVariables, NoMonomorphismRestriction, ViewPatterns #-}
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
import Math.Probably.RandIO
import QuantalHelp
import Baysig.Estimate.RTS hiding (repeat)
import qualified Data.Binary as B
import qualified Numeric.LinearAlgebra as L
import Math.Probably.MCMC
import Math.Probably.FoldingStats
import System.IO
import Data.Ord
import System.Posix.Directory
import System.Cmd

import QuantalHelp

import Graphics.Gnewplot.Exec
import Graphics.Gnewplot.Types
import Graphics.Gnewplot.Style
import Graphics.Gnewplot.Panels
import Graphics.Gnewplot.Instances
import Graphics.Gnewplot.Histogram

import Control.Monad.State.Strict

import Control.Monad.Trans

writeTex fnm doit = do
  h <- openFile fnm WriteMode 
  hPutStrLn h $ unlines     ["\\documentclass[11pt]{article}",
     "%include lhs2TeX.fmt",
     "%include polycode.fmt",
     "\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm]{geometry}",
     "\\usepackage{graphicx}",
     "\\begin{document}"]
  execStateT doit h
  hPutStrLn h "\\end{document}"
  hClose h
  system $ "pdflatex "++fnm
  return ()

putLn s = do
  h <- get 
  lift $ hPutStrLn h s

plotIt nm obj = do
  h <- get 
  lift $ gnuplotToPS (nm++".eps") $ obj
  lift $ system $ "epstopdf "++nm++".eps"
  putLn $"\\includegraphics[width=16cm]{"++nm++"}\n"

simsum = writeTex "simSummary.tex" $ forM_ [1000, 2500] $ \ntr -> do
   let ntrs = pad $ reverse $ drop 3 $ reverse $ show ntr
   pts <- fmap (concat . concat) $ forM [25::Double, 50, 100, 200, 300] $ \ns -> do
            forM [1..5] $ \run -> do
               let nsstr = take 2 $ show $ round ns           
                   sessnm = ntrs++nsstr++show (round run)
                   fnm = sessnm++"/npq_samples"
                   qsim = simq * (100 / realToFrac ns)
               lift $ print fnm
               ifM (lift $ doesFileExist fnm) 
                   (do vsamples::[L.Vector Double] <- lift $ fmap (read)  $ readFile fnm
                       return $ zip (repeat (ns+(run-3)*3)) $ map (\v-> exp (v@>3)) $ thin 100 vsamples)
                   (return [])
   plotIt ("sim_npq_res"++ntrs) pts


whenM mb mu = do 
   b <- mb
   when b mu

ifM mb mp ma = do
   b <- mb
   if b then mp else ma

countSigs sess = do
  putStr $ sess ++ " "
  whenM (doesFileExist (take 6 sess++"/epsps")) $ do
    lns <- fmap lines $ readFile (take 6 sess++"/epsps")
    putStrLn $ show (length lns)
  whenM (doesFileExist (take 6 sess++"/npq_samples")) $ do
    fileconts::[L.Vector Double] <- fmap read $ readFile (take 6 sess++"/npq_samples")
    
    let (mnpars, sdpars) = runStat meanSDF fileconts
    putStrLn $ show $ L.toList mnpars -- ::L.Vector Double)
    putStr $ show $ L.toList sdpars -- ::L.Vector Double)
    let q = exp $ mnpars L.@> 3
    sigs <- do LoadSignals sigs <-  B.decodeFile $ take 6 sess++"/sigs_"++take 6 sess++"_epsps" 
               return sigs
    let wf@(Signal _ _ sv) = baselineSig 0.003 $ averageSigs $ sigs
    let wfAmp = foldl1' max $ L.toList sv
    let qmv  = wfAmp * q
    putStrLn $ "\nq = "++show qmv ++" mV"

    return ()
  putStr "\n"

compareNPQ = writeTex "npqSummary.tex" $ do
  ptss <- forM (zip [(0::Double)..] datasess) $ \(nsess, sess) -> do 
     let fnm = sess++"/npq_samples"
     lift $ print fnm
     ifM (lift $ doesFileExist fnm) 
         (do vsamples::[L.Vector Double] <- lift $ fmap (thin 1000. read)  $ readFile fnm
             sigs <- do LoadSignals sigs <-  lift $ B.decodeFile $ take 6 sess++"/sigs_"++take 6 sess++"_epsps" 
                        return sigs
             let wf@(Signal _ _ sv) = baselineSig 0.003 $ averageSigs $ sigs
             let wfAmp = foldl1' max $ L.toList sv
             let wf1@(Signal _ _ sv1) = baselineSig 0.003 $ averageSigs $ take 20 sigs
             let wfAmp1 = foldl1' max $ L.toList sv1
             return $ (zip (repeat nsess) $ vsamples, (wfAmp, wfAmp1)))
         (return ([], (1, 0)))
  plotIt "nplot" $ ManySup $ map (map (\(i,pts)-> (i, pts L.@> 0))) $ map fst ptss
  plotIt "phiplot" $ ManySup $ map (map (\(i,pts)-> (i, pts L.@> 2))) $ map fst ptss
  let q (iptss, (wfamp, _)) = map (\(i,pts)-> (i, (exp $ pts L.@> 3) * wfamp)) iptss
  plotIt "qplot" $ ManySup $ map q ptss
  let q1 (iptss, (wfamp, amp1)) = map (\(i,pts)-> (i, amp1)) $ take 1 iptss
  plotIt "ampplot" $ ManySup $ map q1 ptss

main = do
  --simsum
--  mapM_ countSigs datasess
 -- compareNPQ
  --plotSigs
  cellSummary "s22mar"

sigZero (Signal dt t0 arr) = Signal dt 0 arr
sigRestrict tstart tend (Signal dt t0 arr) = Signal dt (tstart) $ L.subVector ndrop ntake arr where
  ntake = round $ (tend-tstart)/dt
  ndrop = round $ tstart/dt


plotSigs = writeTex "npqSummary.tex" $ do
 let fls = map (("22Mar"++) . padzero . show) [2..18] --words "22Mar02 22Mar13 22Mar18"
 let sess = "s22mar"
 forM_ fls $ \fl -> do
    LoadSignals sigs <- lift $ B.decodeFile (sess++"/sigs_"++fl++"_epsps")
    putLn fl
    plotIt fl $ map (sigZero) sigs


countSigs1 sess = do
  spks <- fmap (concat . catMaybes) $ inEverySession $ whenContinues sess $ do
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
     return $ Just spikeg
  print $ (sess, length spks)

ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)


cellSummary sess = writeTex "cellSummary.tex" $ do
  (t0s'::[Double], amps'::[Double],sds::[Double]) <- fmap ffile  $ lift $ readFile (sess++"/epsps3")

  (wf, wfAmp, sigs) <- lift $ getWf sess
  let tsamps = filter (getFilter sess) $ zip t0s' $ amps'
      t0s = map fst tsamps


  let weighCurve' = map (weighRegression tsamps ) t0s
      maxPcurve = foldl1 max weighCurve'
      pcurve = map (/(maxPcurve)) weighCurve'
  let globalSd = runStat  meanF sds

  plotIt "pcurve" $ zip t0s pcurve :+: tsamps

  plotIt "wf" wf

  plotIt "wfs" $ concat [take 5 sigs, take 5 $ reverse sigs]
