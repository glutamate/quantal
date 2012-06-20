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

main = do --getArgs >>= dispatch . head
  gen
  main1
dispatch "gen" = gen
dispatch "pdf" = main1


prelease = 0.75
nrel = 25
qreal = 0.01
fnm = "/epsps_"++show nrel++"_"++show (round $ prelease*100)

gen = runRIO $ do
   let covM = mkCovM (np+1) hatPars
       covM1 = mkCovM (np) hatPars
       cholm = chol $ covM
       sess = "calsd"
       tc = 170
       wft t = (((step (t-simt0))*tc*tc*(t-simt0))*(exp ((0.000-(t-simt0))*tc)))
       meanVec = (fillV (np))$(\i-> wft (toD i))
       invDetails = invlndet covM1
   h<- io $ openFile (take 6 sess++fnm) WriteMode 
   forM_ [1..2000] $ \i-> do
        amp <- sample $ binGaussFull nrel prelease qreal 0.2 0
        sig <- sample $ gpByChol dt 
                          (\t-> amp*wft t)
                          cholm
        
        let initialV = L.fromList [0.01,0.1]
        case laplaceApprox defaultAM {nmTol = 0.01} (posteriorSigV (Signal dt 0 meanVec) invDetails sig) [] [] initialV of

          (v, Just cor, smplx) -> do
                let ampHat = v @> 1
                    sd = sqrt $ (L.@@>) cor (1,1)                

                io $ hPutStrLn h $ show (i, ampHat,sd)
                io $ print (i,amp, ampHat, sd)
   io $ hClose h 
   io $ writeFile (take 6 sess++"/noisePars") $ show hatPars
 
post amps v = sum $ (flip map) amps $ \amp->
               binGaussLogPdf (nrel) (prelease) (q) (0.01) (sd) amp
  where sd = v@>0
        q = v@> 1


main1 = do
  let sess = "calsd"
{-  h <- openFile ("FigCalSD.tex") WriteMode 
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
     "\n"] -}
  let ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)

  (t0s'::[Double], amps'::[Double],sds::[Double]) <- fmap ffile  $ readFile (take 6 sess++fnm)

  --plotIt "dalahist" $ Histo 100 amps'
  let glosd = runStat meanF sds
  print glosd

  let lapars@(mid,_,_) = laplaceApprox defaultAM {nmTol = 0.001} (post amps') [] [] $ L.fromList [glosd, 0.01]
  print lapars
  putStr "ratio="
  print $ (mid@>0)/glosd

  putStr "ratioVar="
  print $ ((mid@>0)/glosd)^2

  putStr "npq="
  print $ nrel*prelease*0.01


--  puts "\\end{document}"
  --hClose h

  --system $ "pdflatex FigCalSD.tex"
  return ()

plotSamPdf nm h sam pdf = do
  amps <- sampleNIO 200000 $ sam
  let lo = foldl1' min amps
      hi = foldl1' max amps
      dx = (hi-lo)/400
      scale = (hi-lo)/200 -- not quite right bec histo goes less than spread o
      xs = [lo,lo+dx..hi]
  let plotIt nm obj = do gnuplotToPS (nm++".eps") $ obj
                         system $ "epstopdf "++nm++".eps"
                         hPutStrLn h $"\\includegraphics[width=16cm]{"++nm++"}\n\n"
  
  plotIt nm $ Histo 100 amps :+: (zip xs $ map ((*scale) . pdf) xs)


{-binGauss ns p q cv bgSd = do
     nr  <- binomial ns p
     gaussD (realToFrac nr * q) (sqrt $ q*cv*q*cv*realToFrac nr + bgSd*bgSd) -}

binGaussFull ns p q cv bgSd = do
     nr  <- binomial ns p
     siteAmps <- forM [0..(nr-1)] $ const $ gaussD q (q*cv)
     gaussD (sum siteAmps) (bgSd)

{-binGaussLogPdf ns p q cv bgSd v 
   = log $ bigSum 0 ns $ \nr -> exp $ (normalLogPdf (realToFrac nr * q) (varToTau $ q*cv*q*cv*realToFrac nr + bgSd*bgSd) v) +  (binomialLogProb ns p nr) -}

