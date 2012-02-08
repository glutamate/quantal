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

main = do
  let sess = "57246a"
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

  plotIt "nsigs" $ take 10 sigs

  --plotIt "autoCorr" $ map autoCorrSig $ sigs 

  --plotIt "autoCorrAv" $ [avSigs $ map autoCorrSig $ sigs] 

  vsamples::[L.Vector Double] <- fmap (read) $ readFile (take 6 sess++"/noise_samples")

  plotIt "theta" $ ("theta",  zip [(0::Double)..] $map (exp . (@>0)) vsamples)
  plotIt "sigma" $ ("sigma", zip [(0::Double)..] $map (@>1) vsamples)
  plotIt "obs" $ ("obs", zip [(0::Double)..] $ map  (exp . (@>2)) vsamples)

  plotIt "thetahist" $ ("theta", Histo 50 $ map (exp . (@>0)) vsamples)
  plotIt "sigmahist" $ ("sigma", Histo 50 $ map (@>1) vsamples)
  plotIt "obshist" $ ("obs", Histo 50 $ map  (exp . (@>2)) vsamples)

  puts "posterior predicted noise signals and autocorrelation\n\n"

  fakesigs<-runRIO $ do
    vsample <- sample $ oneOf vsamples
    let sigma = vsample@>1
        theta = exp $ vsample@>0
        obs = exp $ vsample @> 2
    let cholm = chol $ fillM (np+1,np+1) $ \(i,j)-> covOU theta sigma (toD i) (toD j)+ifObs i j obs
    sample $ sequence $ replicate 10 $ gpByChol dt tmax (\t-> 0) cholm
    
  plotIt "fakesigs" $ fakesigs

  fakeautocorr<-runRIO $ sample $ sequence $ replicate 10 $ do
    vsample <- oneOf vsamples
    let sigma = vsample@>1
        theta = exp $ vsample@>0
        obs = exp $ vsample @> 2
    let cholm = chol $ fillM (np+1,np+1) $ \(i,j)-> covOU theta sigma (toD i) (toD j)+ifObs i j obs
    sigs <- sequence $ replicate 20 $ gpByChol dt tmax (\t-> 0) cholm
    return $ avSigs $ map autoCorrSig $ sigs


  plotIt "fakeautoCorr" $ ( fakeautocorr ) :+: Lines [LineWidth 5] [avSigs $ map autoCorrSig $ sigs] 


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