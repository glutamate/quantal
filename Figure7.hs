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

seqit x = x `seq` x
main = do
  h <- openFile ("Figure7.tex") WriteMode 
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

  puts $ "Figure 7"

  hists <- forM ( datasess ) $ \sess -> do
      (wf, wfAmp, _) <- getWf sess
      vsamples::[L.Vector Double] <- fmap (drop (getBurnIn "sess") . map L.fromList . catMaybes . map safeRead . lines) $ readFile (take 6 sess++"/npq_samples")
      let ns =  HistoStyle "histeps" 50 $ map ( (@>0)) vsamples
          ps =  HistoStyle "histeps" 50 $ map (\v->1/(1+exp(-(v @> 2)))) vsamples
          qs =  HistoLog "histeps" 50 $ map ((*wfAmp) . exp . (@>3)) vsamples
          nps = HistoStyle "histeps" 50 $ map (\v -> (v@>2) * (roundD $ (v@>0))) vsamples
      print sess
      print wfAmp
      print (showHisto ns, showHisto ps, showHisto qs)
      return $ (sess, [ns,ps,qs, nps])

  let getN n  = map $ \(sess, vals) -> (take 3 sess, vals!!n)
  let getN' n  = map $ \(sess, vals) -> vals!!n

  --let nsol = 1 -- not really!

  let plotIxN lo hi =     (YLabel "Sol"  $ XRange lo hi $  NoYaxis $ ManySup $ take nsol $ getN 0 hists) 
                       :==: (YLabel "Greg" $ XLabel "N (release sites)" $ XRange lo hi $ NoYaxis $ ManySup $ drop nsol $ getN 0 hists)

  let plotIxP lo hi =     (XRange lo hi $  NoYaxis $ ManySup $ take nsol $ getN' 1 hists) 
                       :==: (XRange lo hi $ XLabel "P" $NoYaxis $ManySup $ drop nsol $ getN' 1 hists)

  let plotIxQ lo hi =     (  LogScaleX $ XRange lo hi $  NoYaxis $ ManySup $ take nsol $ getN' 2 hists) 
                       :==: (  LogScaleX $ XRange lo hi $  XLabel "Q (mV)" $NoYaxis $ ManySup $ drop nsol $ getN' 2 hists)
  let plotIxNP lo hi =     (  ManySup $ take nsol $ getN' 3 hists) 
                        :==: (   ManySup $ drop nsol $ getN' 3 hists)

{-  plotIt "nhists" $ plotIxN 0 800
  plotIt "phists" $ plotIxP 0 1
  plotIt "qhists" $ plotIxQ (-4) (-0.5)
  plotIt "nphists" $ plotIxNP 0 0.07 -}

  plotIt "fig7" $ {- 34 % plotIxN 0 800 :|: 66% (plotIxP 0 1 :||: -}  plotIxQ (1e-4) (0.1)

{-  plotIt "npqhists" $ Hplots [GnuplotBox $ plotIxN 0 800, 
                              GnuplotBox $ plotIxP 0 1, 
                              GnuplotBox $ plotIxQ 0 0.07] -}


  --plotIt "test" $ GnuplotTest

  

  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure7.tex"
  return ()

showHisto (Histo n dbls) = "" -- show (n, dbls)
showHisto (HistoStyle sty n dbls) = "" --show (n, sty, dbls)
showHisto (HistoLog sty n dbls) = "" --show (n, sty, dbls)

safeRead x = case readsPrec 5 x of 
               [] -> Nothing
               (x,_):_ -> Just x



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