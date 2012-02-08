module MPFA where

import Math.Probably.FoldingStats
--import Math.Probably.Distribution
import Math.Probably.GlobalRandoms
import Math.Probably.Sampler
import Math.Probably.StochFun
import Math.Probably.MCMC
import QueryPlots
import QueryTypes
import PlotGnuplot
import TNUtils

import qualified Math.Probably.PDF as PDF
import Control.Monad
import System.IO.Unsafe
epsp :: Double -> Sampler Double
epsp p = do 
  let n =100::Int
      q = 10::Double          
      noise = 1::Double
  nr <- realToFrac `fmap` binomial n p
  gauss (nr*q) (sqrt $ q*nr*noise*noise)
  --return $ q * (realToFrac nr)
  --(q*) `fmap` gauss (n*p) (sqrt $ n*p*(1-p))

--pdfAndSampler :: PDF.PDF Double -> Sampler Double -> GnuplotBox
pdfAndSampler pdf sam = 
    let pts = sampleN 10000 sam
        (min, max) = runStat (both minF maxF) pts
    in FunSeg min max pdf :+: Histo 100 (pts)


tsegs :: Double -> [((Double,Double), ())]
tsegs t = let ts = [0,t..]
          in zip (zip ts (tail ts)) (repeat ())

mpfa :: Double ->   -- segment length
        [(Double,Double)] -> -- events with amplitudes
        [(Double,Double)] -- mean variance points

mpfa t evs = 
  let maxt = maximum $ map fst evs
      segs = map (:[]) $ takeWhile (not . (>maxt) . fst . fst) $ tsegs t
      f seg = case during seg evs of
         [] -> []
         tamps -> [runStat (meanF `both` varF) $ map snd tamps]
  in concatMap f segs 



like1 :: (Double, Double,Double) -> Double -> Double
like1 (n, p, q) i = bin' n p (i/q)
    where bin n p k = (PDF.gaussD (n*p) (sqrt $ n*p*(1-p))) (k)
          bin' n p k = (PDF.binomial (round $ n) p (round $ k))

like2 (n,p,q) i= (PDF.binomial (round $ n) p (round $ i/q))

--((PDF.gaussD (n*p) (sqrt $ n*p*(1-p))) (i/q)) /3



simps, initps :: [Double]
simps = [0.2,0.5,0.8]
initps = map (const 0.5) simps

timeC :: [Double]
timeC =[0,1..1000] 

initparsC, simparsC :: (Double, Double, Double,Double, Int,Double,Double)
simparsC = (0.99, 0.01, 0.02, 500, 100, 10, 1)
initparsC = (0.8, 0.02, 0.02, 500, 40, 25,1 )
initparsC2 = (0.99, 0.01, 0.02, 500, 90, 12, 1)
psC :: [Double]
psC = let pars@(phi, plo, slop, offset, n, q, noise) = simparsC
          ps = map (calcp pars) timeC
       in ps

calcp (phi, plo, slop, offset, n, q, noise) t = (phi-plo)/(1+exp(offset*slop-slop*t))+plo

epspsC :: [Double]
epspsC = head $ unsafePerformIO $ runSamplerIO $ mapM epsp psC

likeC :: PDF.PDF (Double, Double, Double,Double, Int,Double,Double )
likeC pars = 
    sum $ map (like1C pars)  $ zip timeC epspsC

like1C :: (Double, Double, Double,Double, Int,Double,Double ) -> PDF.PDF (Double, Double)
like1C pars@(phi, plo, slop, offset, n, q, noise) (t,v) = 
    let p = calcp pars t
        nrs = [1..n]      
    in log $ sum $ for nrs $ \nr ->  
        let nrdbl = realToFrac nr 
        in PDF.gaussD (nrdbl*q) (sqrt $ q*nrdbl*noise*noise) (v) *
           PDF.binomial n p nr -- (PDF.logBinomial (round $ n) p (round $ v/q))

--for = flip map

proposalC :: (Double, Double, Double,Double, Int,Double,Double ) -> 
             Sampler (Double, Double, Double,Double, Int,Double,Double )
proposalC (phi, plo, slop, offset, n, q, noise) = do
  nphi <- psam' phi
  nplo <- psam' plo
  nslop <- gauss slop 0.00005
  noffset <- gauss offset 1
  nn <- ifProb3 0.2 (n-1) (n+1) n
  nq <- gauss q 0.5
  nnoise <- gauss noise 0.1
  return (phi, plo, slop, offset, nn, nq, noise)

simvs, initvs :: (Double, [Double], Double)
simvs = (100,simps, 10)
initvs = (70, initps, 15)

dat :: [[Double]]
dat = map (sampleN 100 . epsp) simps
dvars = map (runStat varF) dat
dmeans = map (runStat meanF) dat

likeMPFA :: PDF.PDF (Double, [Double], Double)
likeMPFA (n,ps,q) = 
    let predmean = map (*(n*q)) ps
        predvar = map (\p -> n*p*(1-p)*q*q) ps
        dist1 (x1,y1) (x2,y2) = sqrt $ (x2-x1)^^2 + (y2-y1)^^2
        dist = sum $ map (uncurry dist1) $ zip (zip dmeans dvars) (zip predmean predvar)
    in negate $ sum $ map (\x-> x*x) $ map (\(mu,var)-> fitfun (n,q) mu - var) $ zip dmeans dvars

fitfun (n, q) npq = npq*(1-npq/(n*q))*q 
    
pmeans (n,ps,q) =  map (*(n*q)) ps
pvars (n,ps,q) = map (\p -> n*p*(1-p)*q*q) ps

likeData :: PDF.PDF (Double, [Double], Double)
likeData (n,ps,q) = sum $ map (\(pts, p) -> sum $ map (like1 (n,p,q)) pts) $ zip dat ps

proposal (n,ps,q) = do
  nps <- mapM psam' ps
  nn <- gaussD n 1 --ifProb3 0.1 n (n+1) (n-1)
  nq <- gaussD q 0.1 --ifProb3 0.05 q (q+1) (q-1)
  return (nn,nps,nq)

maxLike :: (a -> Double) -> (a -> Sampler a) -> a -> Markov (a, Double)
maxLike like prop init = Mrkv (condSampler sf) init' id
    where sf (x,likeval) = do 
            newx <- prop x
            let newlike = like newx
            return $ if likeval > newlike || nanOrInf newlike
                        then (x,likeval)
                        else (newx, newlike)
          init' = (init, like init)
    
ml = head $ drop 100 $ unsafePerformIO $ runMarkovIO $ maxLike likeData proposal initvs
ml1 = head $ drop 100 $ unsafePerformIO $ runMarkovIO $ maxLike likeData proposal simvs
ml2 = head $ drop 50 $ unsafePerformIO $ runMarkovIO $ maxLike likeMPFA proposal initvs
ml3 = head $ drop 15 $ unsafePerformIO $ runMarkovIO $ maxLike (likeC) proposalC initparsC

--main = putStrLn $ accushow ml4

trpars = \(phi, plo, slop, offset, n, q, noise) -> ((phi, plo), (slop, offset), (n,q, noise))

{-main = do
  print $ likeC simparsC
  mapM (print) =<< (fmap (thin 99 . zip [(0::Int) ..] . map (onFst (accushow . trpars)) ) $ runMarkovIO $ maxLike (likeC) proposalC initparsC2) -}

logit x = 1/(1+exp(negate $ x))
logitInv y = negate $ log (1/y-1)

psam' p = fmap logit $ gauss (logitInv p) 0.05

alpha = 20
newbeta p = round $ (alpha-p*alpha)/p
newbeta' p = (alpha-p*alpha)/p
psam :: Double -> Sampler Double
psam p = do
  beta (round alpha) (newbeta p)

proposalPDF (n,ps, q) (nn,nps,nq) 
    = PDF.gaussD n 1 nn +
      PDF.gaussD q 0.1 nq +
      sum (map (\(p, np) -> (PDF.beta 20 (newbeta' p) np)) (zip ps nps) )


{-beta :: Int -> Int -> Sampler Double
beta a b = 
    let gam n = do us <- forM [1..n] $ const unitSample
                   return $ log $ product us
    in do gama1 <- gam a
          gama2 <- gam a
          gamb <- gam b
          return $ gama1/(gama1+gamb)-}


ifProb :: Double -> a -> a -> Sampler a
ifProb p c a = do
  u <- unitSample
  if u<p then return c
         else return a
ifProb3 :: Double -> a -> a -> a-> Sampler a
ifProb3 p c a1 a2 = do
  u <- unitSample
  case u of
    _ | u < p -> return a1 
    _ | u > (1-p) -> return a2 
    _ | otherwise -> return c 
