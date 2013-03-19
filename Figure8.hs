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
  h <- openFile ("Figure8.tex") WriteMode 
  let puts s = hPutStrLn  h $ s ++ "\n"
      plotIt nm obj = do gnuplotToPS (nm++".eps") $ obj
                         system $ "epstopdf "++nm++".eps"
                         puts $"\\includegraphics[width=16cm]{"++nm++"}\n\n"
      putLn s = 
            lift $ hPutStrLn h s
 
  puts $ unlines     ["\\documentclass[11pt]{article}",
     "%include lhs2TeX.fmt",
     "%include polycode.fmt",
     "\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm]{geometry}",
     "\\usepackage{graphicx}",
     "\\begin{document}",
     "\n"]

  (wf0, _, sigs0) <-  getWf "s22_15"
  plotIt ("wf15") wf0
 
  plotIt ("wfs15") $ take 10 sigs0

  let ffile = (unzip3 .  sortBy (comparing fst3) . map read . lines)

  let wfEpsp@(Signal _ _ sv) = sigFrom 120.83 wf0

  let trueWfAmp = foldl1' max $ L.toList sv

  
  puts $ "Figure 8"
  puts $ "WfAmp = "++show trueWfAmp
  forM (words "s22_05 s22_15 s22_25") $ \sess -> do
    puts $ "cv=0."++ drop 4 sess
    fileconts::[L.Vector Double] <- fmap (map L.fromList . catMaybes . map safeRead . lines) $ readFile (sess++"/npq_samples")
    (wf, wfAmp, sigs) <-  getWf sess




    let ampsFromFile npars = do
        (t0s'::[Double], amps::[Double],sds::[Double]) <- fmap ffile  $ readFile (sess++"/epsps"++show npars)
        let tsamps = zip t0s' $ map (*wfAmp) amps
        return tsamps

    tsamps3 <- ampsFromFile 3

    --plotIt ("tsamps"++sess) tsamps3


    print wfAmp
   
    print $ head fileconts
    print $ last fileconts

    let putIt s f = 
          putLn $ s++" = "++(show $ runStat meanSDF $ map f fileconts)

    let getN = (\v-> (realToFrac $ round $ v @> 0)::Double) 
        getCV= (\v-> exp $ v @> 1)
        getPhi = (\v->1/(1+exp(-(v @> 2))))
        getQ =  (\v->trueWfAmp * (exp $ v @> 3))
    
    --plotIt ("phi"++sess) $ ("phi",zip [0::Double .. ] $ map getPhi fileconts)
    --plotIt ("n"++sess) $ ("n", zip [0::Double .. ] $ map getN fileconts)
    --plotIt ("q"++sess) $ ("q", zip [0::Double .. ] $ map getQ fileconts)


    plotIt ("nhist"++sess) $ ("n", XRange 800 1100 $ Histo 50 $ map getN fileconts)
    plotIt ("qhist"++sess) $ ("q", XRange 0.013 0.016$ Histo 50 $ map getQ fileconts)
    plotIt ("phist"++sess) $ ("pmax", XRange 0.11 0.16 $ Histo 50 $ map getPhi fileconts)


  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure8.tex"
  return ()


safeRead :: Read a => String -> Maybe a
safeRead x = case readsPrec 5 x of 
               [] -> Nothing
               (x,_):_ -> Just x


