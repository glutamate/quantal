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

  vsamples::[L.Vector Double] <- fmap (read) $ readFile (take 6 sess++"/noise_samples3")

  --plotIt "theta" $ ("theta",  zip [(0::Double)..] $map (exp . (@>0)) vsamples)
  --plotIt "sigma" $ ("sigma", zip [(0::Double)..] $map (@>1) vsamples)
  --plotIt "obs" $ ("obs", zip [(0::Double)..] $ map  (exp . (@>2)) vsamples)

  let thetahist =   Histo 50 $ map ( (@>1)) vsamples
  let sigmahist =   Histo 50 $ map (@>0) vsamples
  let obshist = Histo 50 $ map  ( (@>2)) vsamples
  let flathist = Histo 50 $ map  ( (@>3)) vsamples

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

  let mkFakeSigs npars fnm =  do
             vsamples' <- lift $ fmap ( thin 10 .read) $ readFile fnm
             vsample <- fmap ( L.toList ) $ sample $ oneOf vsamples'
             let vsample1 = if npars == 4 then take 2 vsample ++ [(vsample!!2)-2, vsample!!3] else vsample
             let cholm = L.chol $ mkCovM (np+1) $ take npars $ vsample1
             sample $ sequence $ replicate 10 $ fmap (id {-baselineSig 0.1 -}) $ gpByChol dt (\t-> 0) cholm
    

  let mkFakeAutoCorr npars fnm = runRIO $ sequence $ replicate 40 $ do
                         sigs <-  mkFakeSigs npars fnm                          
                         return (avSigs $ map autoCorrSig $ sigs, map sigSpan sigs) 



  let mkAutoCorrPlot fake = ( Lines [LineWidth 1, LineType 1, LineColor "red"] $ map fst fake ) :+: Lines [LineWidth 5, LineType 1, LineColor "black"] [avSigs $ map autoCorrSig $ sigs] 


  fakesigs4 <- runRIO $ mkFakeSigs 4 (take 6 sess++"/noise_samples4")
  fakeautocorr4 <- mkFakeAutoCorr 4 (take 6 sess++"/noise_samples4")

  fakesigs3 <- runRIO $  mkFakeSigs 3 (take 6 sess++"/noise_samples3")
  fakeautocorr3 <- mkFakeAutoCorr 3 (take 6 sess++"/noise_samples3")

  fakesigs2 <- runRIO $  mkFakeSigs 2 (take 6 sess++"/noise_samples2")
  fakeautocorr2 <- mkFakeAutoCorr 2 (take 6 sess++"/noise_samples2")

  fakesigs1 <- runRIO $  mkFakeSigs 1 (take 6 sess++"/noise_samples1")
  fakeautocorr1 <- mkFakeAutoCorr 1 (take 6 sess++"/noise_samples1")

  --plotIt "plotr" $ realNoise :||: Noplot


  --plotIt "plot4" $ fakesigs4 :||: mkAutoCorrPlot fakeautocorr4
  --plotIt "plot3" $ fakesigs3 :||: mkAutoCorrPlot fakeautocorr3
  --plotIt "plot2" $ fakesigs2  :||: mkAutoCorrPlot fakeautocorr2
  --plotIt "plot1" $ fakesigs1 :||: mkAutoCorrPlot fakeautocorr1



  --plotIt "fig2" $ ((realNoise :==: fakeNoise3) :||: autoCorr3) :==: lowplot
  --plotIt "fig2" $ (plot3v  (C $ mkAutoCorrPlot fakeautocorr3) (Aii fakesigs1) (Ai realNoise)) 
  --                 :||: (B $ plot3v ("theta", thetahist) ("sigma", sigmahist) ("observation", obshist))

  --plotIt "sigma" $ ("sigma", zip [(0::Double)..] $map (@>0) vsamples)
  --plotIt "theta" $ ("theta",  zip [(0::Double)..] $map (@>1) vsamples)
  --plotIt "obs" $ ("obs", zip [(0::Double)..] $ map  (@>2) vsamples)
  --plotIt "flat" $ ("flat", zip [(0::Double)..] $ map  (exp . (@>3)) vsamples)

  puts $ "sigspan real = "++ show (runStat meanF $ map sigSpan realNoise) ++"\n\n"
  --puts $ "sigspan model1 = "++ show (runStat meanSDF $ map (runStat meanF . snd) fakeautocorr1) ++"\n\n"
  --puts $ "sigspan model2 = "++ show (runStat meanSDF $ map (runStat meanF . snd) fakeautocorr2) ++"\n\n"
  puts $ "sigspan model3 = "++ show (runStat meanSDF $ map (runStat meanF . snd) fakeautocorr3) ++"\n\n"
  --puts $ "sigspan model4 = "++ show (runStat meanSDF $ map (runStat meanF . snd) fakeautocorr4) ++"\n\n"

  let spanhists = ("1", Histo 20 (map (runStat meanF . snd) fakeautocorr1)) :+:
                       ("2", Histo 20 (map (runStat meanF . snd) fakeautocorr2)) :+:
                       ("3", Histo 20 (map (runStat meanF . snd) fakeautocorr3)) :+:
                       ("4", Histo 20 (map (runStat meanF . snd) fakeautocorr4)) 


  scatters <- fmap (unzip3) $ forM datasess $ \dsess -> do
     vsams3::[L.Vector Double] <-  readFile' (take 6 dsess++"/noise_samples3")
     vsams2::[L.Vector Double] <-  readFile' (take 6 dsess++"/noise_samples2")
     return (( zip (map (@>0) vsams2) (map (exp . (@>1)) vsams2)), 
             (dsess, zip (map (@>0) vsams3) (map (exp . (@>1)) vsams3)), 
             (zip (map (exp . (@>2)) vsams3) (map (exp . (@>1)) vsams3)))
     
  --plotIt "scatpars2" $ ManySup $ fst3 scatters
  --plotIt "scatpars3" $ ManySup $ snd3 scatters
  --plotIt "scatpars4" $ ManySup $ trd3 scatters
  {- 0: sigma, 1: logtheta, 2: logobs -}

  --let datareal = realNoise
  let triplot a b c = 33% a :|: 67%(b :||: c)


  let topRow = 33% realNoise :|: 67%(fakesigs2 :||: fakesigs3)
  let midRow = triplot (ManySup $ fst3 scatters ) (ManySup $ snd3 scatters) (ManySup $ trd3 scatters)
  let botPlot = spanhists :||: mkAutoCorrPlot fakeautocorr3

  plotIt "finalFigure1" $ (A realNoise :||: B fakesigs1) :==: (C fakesigs2 :||: D fakesigs3)
  plotIt "finalFigure2" $ ((Ei $ YRange 0 500 $ XRange 1 7 $ ManySup $ fst3 scatters ) :||: 
                           (Eii $ YRange 0 500 $ XRange 1 7 $ ManySup $ snd3 scatters)) :==: 
                          ((Eiii $ YRange 0 500 $ ManySup $ trd3 scatters) :||: 
                          (Fi spanhists :==: Fii (mkAutoCorrPlot fakeautocorr3)))

  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure2.tex"
  return ()



readFile' fnm = do
  ex <- doesFileExist fnm
  if ex 
     then fmap ( thin 10 .read) $ readFile fnm
     else return []

sigSpan :: Signal Double -> Double
sigSpan (Signal dt t0 vec) = let (lo,hi) = runStat (both minF maxF) $ L.toList vec
                             in hi - lo

autoCorrSig :: Signal Double -> Signal Double
autoCorrSig (Signal dt t0 vec) = Signal dt 0 $ L.fromList $ autoCorr $ L.toList vec

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
