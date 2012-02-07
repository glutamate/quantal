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
  let sess = "00c9bd"
  h <- openFile ("Figure4.tex") WriteMode 
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

  amps <- sampleNIO 50000 $ binGaussFull 10 0.2 0.1 0.2 0.01
  let lo = foldl1' min amps
      hi = foldl1' max amps
      dx = (hi-lo)/400
      scale = (hi-lo)/200 -- not quite right bec histo goes less than spread o
      xs = [lo,lo+dx..hi]

  let pdf = binGaussPdf 10 0.2 0.1 0.2 0.01

  plotIt "hist1" $ Histo 100 amps :+: (zip xs $ map ((*scale) . pdf) xs)

  puts "\\end{document}"
  hClose h

  system $ "pdflatex Figure4.tex"
  return ()

plotSamPdf h sam pdf = undefined

binGauss ns p q cv bgSd = do
     nr  <- binomial ns p
     gaussD (realToFrac nr * q) (sqrt $ q*cv*q*cv*realToFrac nr + bgSd*bgSd)

binGaussFull ns p q cv bgSd = do
     nr  <- binomial ns p
     siteAmps <- forM [0..(nr-1)] $ const $ gaussD q (q*cv)
     gaussD (sum siteAmps) (bgSd)

{-binGaussLogPdf ns p q cv bgSd v 
   = log $ bigSum 0 ns $ \nr -> exp $ (normalLogPdf (realToFrac nr * q) (varToTau $ q*cv*q*cv*realToFrac nr + bgSd*bgSd) v) +  (binomialLogProb ns p nr) -}

