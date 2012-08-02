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
import qualified Math.Probably.PDF as PDF
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
import qualified Data.Array.IArray as IA


main = do
 let sess = "s22mar"
 let fls = words "22Mar02 22Mar13 22Mar18"
 let t0s = map (*60) [0, 18,38]
 forM_ (zip fls t0s) $ \(fl,t0) -> do
   lns <- fmap lines $ readFile $ "/home/tomn/Dropbox/Quantal/Simmons/"++fl++".txt"
   mapM print $ take 40 lns
   let (vm, more) = (unParser getEpsp) lns
   let (map (+t0) . init . tail -> trigs, _) = (unParser getTriggers) more
 --mapM print $ take 10 trigs
   print $ (realToFrac (length vm) * 0.0002)
   let vmSig = Signal 0.0002 t0 $ L.fromList vm
-- encodeFile (fl++"/bigvm") $ LoadSignals [vmSig]
   --print "done sigs"
   encodeFile (sess++"/sigs_"++fl++"_epsps") $ LoadSignals (map (\t->limitSig (t-0.03) (t+0.17) vmSig) $ trigs)
   encodeFile (sess++"/sigs_"++fl++"_noise") $ LoadSignals (map (\t->limitSig (t-0.23) (t-0.03) vmSig) $ tail trigs)
 writeFile (sess++"/sessions") $ show fls
 return ()

getEpsp :: Parser String [Double]
getEpsp = do
  dropTillP ("\"START" `isPrefixOf`) 
  repeatUntil (=="\r") read headP

limitSig lo hi (Signal dt t0 arr) = 
    let ndrop = round $ (lo - t0)/dt
        ntake = round $ (hi - lo)/dt
    in Signal dt lo (L.subVector ndrop ntake arr)


getTriggers :: Parser String [Double]
getTriggers = do
  dropTillP ("\"Trig" `isPrefixOf`) 
  _ <- headP
  repeatUntil (=="\r") read headP

  
newtype Parser s a = Parser { unParser :: [s] -> (a,[s]) }

instance Monad (Parser s) where
   return x = Parser $ \ss -> (x,ss)
   mx >>= f = Parser $ \ss -> let (x, ss') = unParser mx ss
                              in unParser (f x) ss'

headP = Parser $ \(s:ss) -> (s, ss)

takeP n = Parser $ \(ss) -> splitAt n ss
            
dropTillP p = Parser go where
   go [] = ((), [])
   go (s:ss) | p s = ((), ss)
             | otherwise = go ss

repeatUntil :: Monad m => (a -> Bool) -> (a -> b) -> m a -> m [b]
repeatUntil p f mx = do x <- mx
                        if p x 
                           then return []
                           else liftM2 (:) (return (f x)) (repeatUntil p f mx)

