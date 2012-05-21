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

plot3v x y z = Vplots [GnuplotBox x, GnuplotBox y, GnuplotBox z]
plot3h x y z = Hplots [GnuplotBox x, GnuplotBox y, GnuplotBox z]


main = do
  let sess = "00c9bd"
  h <- openFile ("Figure2.tex") WriteMode 
  let puts = hPutStrLn h
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

  LoadSignals sigs <- decodeFile $ take 6 sess++"/sigs_"++take 6 sess++"_noise"

  let realNoise =  take 10 sigs

  --plotIt "autoCorr" $ map autoCorrSig $ sigs 

  --plotIt "autoCorrAv" $ [avSigs $ map autoCorrSig $ sigs] 

  vsamples::[L.Vector Double] <- fmap (drop 1000 . read) $ readFile (take 6 sess++"/noise_samples")

  --plotIt "theta" $ ("theta",  zip [(0::Double)..] $map (exp . (@>0)) vsamples)
  --plotIt "sigma" $ ("sigma", zip [(0::Double)..] $map (@>1) vsamples)
  --plotIt "obs" $ ("obs", zip [(0::Double)..] $ map  (exp . (@>2)) vsamples)

  let thetahist =   Histo 50 $ map (exp . (@>0)) vsamples
  let sigmahist =   Histo 50 $ map (@>1) vsamples
  let obshist = Histo 50 $ map  (exp . (@>2)) vsamples

  puts "posterior predicted noise signals and autocorrelation\n\n"

  {-hists <- forM datasess $ \sess -> do
      (wf, wfAmp, _) <- getWf sess
      vsamples::[L.Vector Double] <- fmap (drop (getBurnIn "sess") . thin 10 . read) $ readFile (take 6 sess++"/noise_samples")
      let ts =  HistoStyle "histeps" 50 $ map (exp . (@>0)) vsamples
          ss =  HistoStyle "histeps" 50 $ map (@>1) vsamples
          os =  HistoStyle "histeps" 50 $ map (exp . (@>2)) vsamples
      print sess
      print wfAmp
      print (showHisto ts, showHisto ss, showHisto os)
      return $ (sess, [ts,ss,os])

  let getN' n  = map $ \(sess, vals) -> vals!!n
  
  let plotIx ix lo hi =     (  ManySup $ take nsol $ getN' ix hists) 
                         :==: (   ManySup $ drop nsol $ getN' ix hists)

  let lowplot = plot3h (plotIx 0 0 1) (plotIx 1 0 1) (plotIx 2 0 1) -}

  fakesigs<-runRIO $ do
    vsample <- sample $ oneOf vsamples
    let sigma = vsample@>1
        theta = exp $ vsample@>0
        obs = exp $ vsample @> 2
        covM = fillM (np+1,np+1) $ \(i,j)-> covOU theta sigma (toD i) (toD j)+ifObs i j obs
    let cholm = L.trans $ L.chol $ covM
    sample $ sequence $ replicate 10 $ gpByChol dt tmax (\t-> 0) cholm
    
--  plotIt "fakesigs" $ fakesigs 

  let fakeNoise3 = map (baselineSig 0.1)  fakesigs
  let fakeNoise2 = Noplot

{-  covM <-runRIO $ do
    vsample <- sample $ oneOf vsamples
    let sigma = vsample@>1
        theta = exp $ vsample@>0
        obs = exp $ vsample @> 2
        covM = fillM (np+1,np+1) $ \(i,j)-> covOU theta sigma (toD i) (toD j)+ifObs i j obs
    return covM

  print $ L.takeDiag covM -}

  fakeautocorr<-runRIO $ sample $ sequence $ replicate 10 $ do
    vsample <- oneOf vsamples
    let sigma = vsample@>1
        theta = exp $ vsample@>0
        obs = exp $ vsample @> 2
    let cholm = chol $ fillM (np+1,np+1) $ \(i,j)-> covOU theta sigma (toD i) (toD j)+ifObs i j obs
    sigs <- sequence $ replicate 20 $ gpByChol dt tmax (\t-> 0) cholm
    return $ avSigs $ map autoCorrSig $ sigs


  let autoCorr3 = ( Lines [LineWidth 1, LineType 1, LineColor "red"]fakeautocorr ) :+: Lines [LineWidth 5, LineType 1, LineColor "black"] [avSigs $ map autoCorrSig $ sigs] 
 

  --plotIt "fig2" $ ((realNoise :==: fakeNoise3) :||: autoCorr3) :==: lowplot
  plotIt "fig2" $ (plot3v  (C autoCorr3) (Aii fakeNoise3) (Ai realNoise)) 
                   :||: (B $ plot3v ("theta", thetahist) ("sigma", sigmahist) ("observation", obshist))


  plotIt "theta" $ ("theta",  zip [(0::Double)..] $map (exp . (@>0)) vsamples)
  plotIt "sigma" $ ("sigma", zip [(0::Double)..] $map (@>1) vsamples)
  plotIt "obs" $ ("obs", zip [(0::Double)..] $ map  (exp . (@>2)) vsamples)


  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure2.tex"
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

showHisto (Histo n dbls) = show (n, dbls)
showHisto (HistoStyle sty n dbls) = show (n, sty, dbls)
