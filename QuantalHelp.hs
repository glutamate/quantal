{-# LANGUAGE DeriveDataTypeable, ScopedTypeVariables, NoMonomorphismRestriction, ViewPatterns, PackageImports, BangPatterns #-}
module QuantalHelp where
import "probably" Math.Probably.MCMC
import "probably" Math.Probably.StochFun
import "probably" Math.Probably.Sampler
import "probably" Math.Probably.FoldingStats
import Control.Monad
--import Data.Array
--import qualified Statistics.Math as SM
import qualified Math.Probably.Student as S
import qualified Data.StorableVector as SV
--import Codec.Image.DevIL
import qualified Data.Array.Unboxed as UA
import Foreign.ForeignPtr
import Foreign.Storable.Tuple
import qualified Data.StorableVector as SV
import System.IO.Unsafe
import qualified Math.Probably.PDF as PDF
import qualified Numeric.LinearAlgebra as L
import "baysig" Baysig.Estimate.RTS
import "probably" Math.Probably.RandIO
import "probably" Math.Probably.NelderMead
import Data.List
import System.Environment
import Data.Array.IArray
import System.Directory


import Data.Binary

import Control.Spoon

import Query hiding (io) 
import QueryTypes
import QueryUtils hiding (averageSigs)

import Debug.Trace

neg = \x-> 0.000-x
square = \x-> x*x
step = \x-> if (x<0.000) then 0.000 else 1.000
fac = \_arg0-> case _arg0 of 1 -> 1; n -> n*(fac (n-1))
assert_fac3 = (fac 3)==6
assert_pairs0 = 2==(fst ((2,3)))
assert_pairs1 = 3==(fst ((2,3)))
assert_case = case 1 of 0 -> False; 1 -> True
--oneOf = \xs-> ((fmap floor)$(uniform 0.000 (realToFrac$(length xs))))>>=(\idx-> return$(ix idx xs))
uniformLogPdf = \lo-> \hi-> \x-> 
   if ((x<hi)&&(x>lo)) 
       then 1.000 
       else {-trace ("expected "++show x++"to be in rng "++show (lo,hi)) -}(0.000-1.000e30)
uniform' = \lo-> \hi-> unit>>=(\x-> return ((x*(hi-lo))+lo))
unormal' = unit>>=(\u1-> unit>>=(\u2-> return ((sqrt ((0.000-2.000)*(log u1)))*(cos ((2.000*pi)*u2)))))
oneTo = \n-> (uniform 0.500 ((realToFrac n)+0.500))>>=(\x-> return (round x))
oneToLogPdf = \hi-> \x-> if ((x<(hi+1))&&(x>0)) then 1.000 else (0.000-1.000e50)
normal = \mean-> \tau-> unormal>>=(\u-> return ((u/(sqrt (tau*2.000)))+mean))
normalSd = \mean-> \sd-> unormal>>=(\u-> return ((u*sd)+mean))
normalPdf = \mu-> \tau-> \x-> (sqrt ((tau/2.000)*pi))*(exp (0.000-(((x-mu)*(x-mu))*tau)))
normalLogPdf :: Double -> Double -> Double -> Double
normalLogPdf = \mu-> \tau-> \x-> (log (sqrt ((tau/2.000)*pi)))+(0.000-(((x-mu)*(x-mu))*tau))
sdToTau = \sd-> 1.000/((2.000*sd)*sd)
meanCVToTau = \mn-> \cv-> 1.000/((2.000*(cv*mn))*(cv*mn))
varToTau :: Double -> Double
varToTau = \var-> recip $ 2*var
tauToSD = \t-> 1.000/(sqrt (t*2.000))
logNormal = \mu-> \tau-> (fmap exp)$(normal mu tau)
logNormalPdf = \mu-> \tau-> \x-> ((sqrt tau)*(1.000/x))*(exp (0.000-(((tau*((log x)-mu))*((log x)-mu))/2.000)))
binNormalApprox = \n-> \p-> normal (n*p) (varToTau ((n*p)*(1.000-p)))
binomialProb = \n-> \p-> \k-> ((choose n k)*(p^k))*((1.000-p)^(n-k))
binomialLogProb :: Int ->Double ->Int -> Double 
binomialLogProb = \n-> \p-> \k-> ((log$(choose n k))+((realToFrac k)*(log p)))+((realToFrac (n-k))*(log (1.000-p)))

binomialLogProb' n p k | k > n = binomialLogProb' n p n 
                       | k < 0 = binomialLogProb' n p 0
                       | n < 0 = error "n zero"
                       | otherwise = ((log$(choose n k))+((realToFrac k)*(log p)))+((realToFrac (n-k))*(log (1.000-p)))

bigLogExpSum lo hi f =  log_sum_exp $ map f [lo..hi]


log_sum_exp xs =  

  -- Choose c to be the element of v that is largest in absolute value.
  let max = runStat maxF xs
      maxabs = runStat (before maxF abs) xs
      c = max {-if ( maxabs > max ) 
             then runStat minF xs
             else max -}

  in log(sum  $ map (\x-> exp(x-c)) xs) + c


countTrue = \_arg0-> case _arg0 of Nil  -> 0; Cons True  bs -> 1+(countTrue bs); Cons False  bs -> countTrue bs
betaLogPdf = \a-> \b-> \x-> log$(((1.000/(S.beta (realToFrac a) (realToFrac b)))*(x^(a-1)))*((1.000-x)^(b-1)))
unfoldN = \n-> \m-> \lastx-> \s-> if (n<m) then ((s n lastx)>>=(\x-> (((unfoldN (n+1) m) x) s)>>=(\xs-> return (Cons x xs)))) else ((s n lastx)>>=(\v-> return (Cons v Nil)))
unfold = \n-> \lastx-> \s-> ((unfoldN 1 n) lastx) s
binGauss = \ns-> \p-> \q-> \cv-> \bgSd-> (binomial ns p)>>=(\nr-> normal ((realToFrac nr)*q) (varToTau$(((((q*cv)*q)*cv)*(realToFrac nr))+(bgSd*bgSd))))
binGaussPdf = \ns-> \p-> \q-> \cv-> \bgSd-> \v-> (bigSum 0 ns)$(\nr-> ((normalPdf ((realToFrac nr)*q) (varToTau$(((((q*cv)*q)*cv)*(realToFrac nr))+(bgSd*bgSd)))) v)*((binomialProb ns p) nr))

--binGaussPdfFrom1 = \ns-> \p-> \q-> \cv-> \bgSd-> \v-> (bigSum 1 ns)$(\nr-> ((normalPdf ((realToFrac nr)*q) (varToTau$(((((q*cv)*q)*cv)*(realToFrac nr))+(bgSd*bgSd)))) v)*((binomialProb ns p) nr))

-- log(x + y) = log(x) + log(1 + exp(log(y) - log(x)))
-- http://www.perlmonks.org/index.pl?node_id=974222
--binGaussLogPdf' = \ns-> \p-> \q-> \cv-> \bgSd-> \v-> log$((bigSum 0 ns)$(\nr-> exp$(((normalLogPdf ((realToFrac nr)*q) (varToTau$(((((q*cv)*q)*cv)*(realToFrac nr))+(bgSd*bgSd)))) v)+((binomialLogProb ns p) nr))))

binGaussLogPdf :: Int -> Double -> Double -> Double -> Double -> Double -> Double 
binGaussLogPdf ns p q cv  bgSd v = 
   bigLogExpSum 0 ns $ \nr-> 
       normalLogPdf ((realToFrac nr)*q) (varToTau$ f1*(realToFrac nr)+f2) v
     + binomialLogProb ns p nr
 where f1 = q*cv*q*cv
       f2 = bgSd*bgSd
        

--binGaussFrom1LogPdf ns p q  cv bgSd v = log$((bigSum 1 ns)$(\nr-> exp$(((normalLogPdf ((realToFrac nr)*q) (varToTau$(((((q*cv)*q)*cv)*(realToFrac nr))+(bgSd*bgSd)))) v)+((binomialLogProb ns p) nr))))


normalInit = \mu-> \(_) -> mu
uniformInit = \lo-> \hi-> (hi+lo)/2.000
oneToInit = \n-> div n 2
binGaussInit = \n-> \p-> \q-> \(_) -> \(_) -> ((realToFrac n)*p)*q
betaInit = \a-> \b-> a/(a+b)
alpha = \tc-> \t-> ((((step t)*tc)*tc)*t)*(exp ((0.000-t)*tc))
qsig = \amp-> \tc-> \t0-> \off-> \t-> off+(amp*(alpha tc (t-t0)))
covOUOld = \theta-> \sigma-> \s-> \t-> (((sigma*sigma)*0.500)/theta)*(exp (0.000-(theta*(abs (s-t)))))
dt = 5.000e-5 -- 0.0002 -- 5.000e-5
tmax = 0.1 --0.2
np = round$(tmax/dt)
toD = \i-> (realToFrac i)*dt

ifObs = \i-> \j-> \sig-> if (i==j) then sig else 0.000

gpByInvLogPdf = \(_) -> \(_) -> \meansig-> \lndet-> \covinv-> \obssig-> let ((dt,_),obsvec) = observe obssig; meanVec = (fillV np)$(\i-> meansig (toD i)) in ((mvnPdf lndet covinv) meanVec) obsvec

posteriorNoiseV 1 sigs v = 
  let inv = L.diag $ L.buildVector np' $ const $ recip $ exp logvar
      lndet = realToFrac np' *  logvar 
  in dens (inv,lndet)
  where logvar = v@> 0
        Signal _ _ sigv1 : _ = sigs
        np' =  L.dim sigv1
        tmax' = realToFrac np' * dt 
        --means = L.toList $ L.subVector 1 10 v
        vmean = last $ L.toList $ v
        dens (inv,lndet) = 
--                uniformLogPdf (0.000) (10.000) sigma
                
-- +uniformLogPdf (0.000-80.000) (0.000-40.000) vmean
                (sum $ (flip map) (zip [1..10] sigs) $ \(i, sigv)->gpByInvLogPdf (dt) (tmax') (\y-> vmean) (lndet) (inv) sigv)

posteriorNoiseV npars sigs v = 
  let covM= mkCovM np' $ take npars $ L.toList v 
  in case spoon $ invlndetC covM of
                        Just (inv,lndet) ->  dens (inv,lndet)
                        _ -> -1000000 -- error $"invlndet: "++show v
  
  where sigma = v@> 0
        Signal _ _ sigv1 : _ = sigs
        np' =  L.dim sigv1
        tmax' = realToFrac np' * dt 
        vmean = last $ L.toList $ v
        dens (inv,lndet) = 
                uniformLogPdf (0.000) (10.000) sigma
                
-- +uniformLogPdf (0.000-80.000) (0.000-40.000) vmean
                +(sum $ (flip map) (zip [1..10] sigs ) $ \(i, sigv)->gpByInvLogPdf (dt) (tmax') (\y-> vmean) (lndet) (inv) sigv)


sigSub x (Signal dt t0 vec) = Signal dt t0 $ L.mapVector (\val -> val - x) vec 

--vOU :: Double -> DO
covOU sigma logtheta ij = (sigma*sigma*0.500/theta)*(exp (-theta*(abs (realToFrac ij*dt))))
  where theta = exp logtheta


mkCovM np' [sigma,logtheta,logobs]  = 
  let h ij = covOU sigma logtheta ij + if ij==0 then (exp logobs) else 0.000 

      line1V :: Array Int Double = listArray (0, np'-1) $ map h [0..np'-1]

      covm' = fillM (np',np') $ \(i,j)-> line1V ! (abs $ i-j)
  in covm' 

mkCovM np' [sigma,logtheta]  = 
  let h ij =  covOU sigma logtheta ij 

      line1V :: Array Int Double = listArray (0, np'-1) $ map h [0..np'-1]

      covm' = fillM (np',np') $ \(i,j)-> line1V ! (abs $ i-j)
  in covm' 

mkCovM np' [logvar]  = 
  let covm' = fillM (np',np') $ \(i,j)-> if i==j then exp logvar else 0
  in covm' 


mkCovM np' [sigma,logtheta,logobs,smoothsd] = 
  let line1V = L.buildVector np' $ \(ij) -> covOU sigma logtheta ij + if ij==0 then (exp logobs) else 0.000 

      wdiff = L.buildVector 20 $ \(j) -> exp $ PDF.gauss 0 (exp smoothsd) (realToFrac j * dt)
      

      line4V :: Array Int Double = listArray (0, np'-1) $ map h [0..np'-1]

      h i =
        let otherixs = [max 0 (i-4)..min (np'-1) (i+4)] 
            (v, wsum) = foldl' (\(!accv, !accw) (nv, nw) -> (nv*nw+accv, nw+accw)) (0,0) $ map (g i) otherixs

        in v/wsum
  
{-  f i j | j < 0 = (0,0) 
        | j >= np = (0,0)
        | otherwise = (line1V L.@> j, wdiff L.@> (abs $ i-j)) -}

      g i j = (line1V L.@> j, wdiff L.@> (abs $ i-j))


      covm' = {-trace (show (take 10 $ L.toList wdiff) ++ "\n" ++ show (take 10 $ L.toList line1V) ++ "\n" ++ show (take 10 $ elems line4V)) $ -} fillM (np',np') $ \(i,j)-> line4V ! (abs $ i-j)
  in covm' 

unSig (Signal _ _ sigv1) = sigv1

sumVec = L.foldVector (+) 0

(@>) = (L.@>)

--'luSolve' . 'luPacked'

invlndetC :: Matrix Double -> (Matrix Double, Double)
invlndetC m = (im, lndet)
  where c = L.chol  m
        im = L.cholSolve c $ L.ident (L.rows m)
        prod = L.foldVector (\x acc -> acc+ log x) 0 $ L.takeDiag c
        lndet = 2* prod 

--  covm = fillM (np',np') $ \(i,j)-> line1V @> (abs $ i-j)

ouSynapseLogPdf (covinv, lndet) 
      = \meansig-> \obssig-> let (_,obsvec) = observe obssig
                                 (_,meanvec) = observe meansig
                             in ((mvnPdf lndet covinv) meanvec) obsvec

scaleSig off x (Signal dt t0 vec) = Signal dt t0 (L.mapVector (\v-> v*x+off) vec)

sigNpts (Signal _ _ v)= L.dim v

zeroSigInit :: Double -> Signal Double -> Signal Double
zeroSigInit tbase (Signal dt tst vec) = 
   let nzero = round $ tbase/dt
       basevec = L.subVector (nzero) (L.dim vec - nzero) vec
       zerovec = L.fromList $ replicate nzero 0
   in Signal dt tst (L.join [zerovec, basevec])

baselineSig :: Double -> Signal Double -> Signal Double
baselineSig tbase (Signal dt tst vec) = 
   let ntake = round $ tbase/dt
       basevec = L.subVector 0 ntake vec
       xsub = (L.foldVector (+) 0 basevec) / realToFrac ntake      
   in (Signal dt tst (L.mapVector (subtract xsub) vec))

baselineSigFrom :: Double -> Double -> Signal Double -> Signal Double
baselineSigFrom tfrom tbase (Signal dt tst vec) = 
   let ntake = round $ (tbase-tfrom)/dt
       basevec = L.subVector (round $ tfrom/dt) ntake vec
       xsub = (L.foldVector (+) 0 basevec) / realToFrac ntake      
   in (Signal dt tst (L.mapVector (subtract xsub) vec))

cacheIn fnm mx = do
  ex<- doesFileExist fnm
  if ex 
     then fmap read $ readFile fnm
     else do x <- mx 
             writeFile fnm $ show x
             return x


initAdaMetFromCov nsam pdf initv retries cov = do
  io $ putStrLn $ "starting from existing cov; try number "++ show retries
  iniampar <- sample $ initialAdaMetFromCov nsam (pdf) initv
                                                 (PDF.posdefify $ cov) 
  io $ print iniampar
  let rate = realToFrac (count_accept iniampar) / realToFrac nsam
  case () of
     _ | retries > 8 -> do io $ putStrLn "initals ok." 
                           return iniampar
     _ | rate > 0.5 -> initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 2 cov 
     _ | rate > 0.40 -> initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 1.5 cov 
     _ | rate < 0.025 ->  initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 0.1 cov 
     _ | rate < 0.075 ->  initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 0.2 cov 
     _ | rate < 0.15 ->  initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 0.3 cov 
     _ | rate < 0.19 ->  initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 0.5 cov 
     _ | rate < 0.24 ->  initAdaMetFromCov nsam pdf initv (retries +1) $ L.scale 0.8 cov 
     _ | otherwise -> do io $ putStrLn "initals ok." 
                         return iniampar
  

-- 

posteriorSigV wf invDetails sig v 
  = ouSynapseLogPdf invDetails (scaleSig (v0) amp wf) sig
 where v0 = v@> 0
       amp = v@>1       

posteriorNPQV amps pcurve sd v = -- ((n,cv,slope,offset,phi,plo,q,tc,t0), loopvals) = 
 oneToLogPdf (800) n
 + normalLogPdf (-1.9) 100 (log cv) -- uniformLogPdf 0.00001 0.5 cv
 +uniformLogPdf (0) (1) phi
 +uniformLogPdf (0.000) (1) q
 +(sum $ (flip map) (zip pcurve amps) $ \(pcurveVal, amp)->
               let p=pcurveVal * phi 
               in binGaussLogPdf (n) (p) (q) (cv) (sd) amp)
  where n = round $ v @> 0
        cv = exp $ v @> 1
        phi = v @> 2
        q = exp $ v @> 3

type Vec = L.Vector Double

fastNPQ :: (Int -> Vec -> Double) -> Int -> Vec -> (Vec, Int, Double, Simplex)
fastNPQ pdfN n0 par0 = fN initLike (n0-1) par0 where
  initLike = -1e20 -- fst $ optimise n0 par0

  fN lastLike nlast pars = let (thislike, thispars, simp) = optimise (nlast+1) pars
                           in if thislike > lastLike 
                                 then fN thislike (nlast+1) thispars
                                 else (pars, nlast, lastLike, simp)
  optimise n pars = let pdf = pdfN n
                        (maxV, _,smplx) = laplaceApprox defaultAM {nmTol = 5} pdf [] [] pars
                        --(maxPost,hess) = hessianFromSimplex (negate . pdfN) [] $ augmentSimplex n smplx 
                        like = pdf maxV
                      in  trace ("laplace"++show n++": "++show maxV++" lh="++show like) $ (like, maxV, smplx)


sigAv (Signal _ _ lv) = runStat meanF $ L.toList lv
       
posDefCov ampar = ampar { ampCov = PDF.posdefify $ ampCov ampar}

augmentSimplex n = map f where
   f (vec, like) = (L.join [morev, vec], like)
   morev = L.fromList [realToFrac n]
       


                  
hessianFromSimplex' :: (L.Vector Double -> Double) -> [Int] -> Simplex -> (L.Vector Double, Matrix Double)
hessianFromSimplex' f isInt sim = 
  let mat :: [L.Vector Double]
      mat = L.toRows $ L.fromColumns $ map fst sim
      fsw ((y0, ymin),ymax) = (y0, max (ymax-y0) (y0-ymin))
      swings = flip map mat $ runStat (fmap fsw $ meanF `both` minFrom 1e80 `both` maxFrom (-1e80)) . L.toList 
      n = length swings
      xv = L.fromList $ map fst swings
      fxv = f  xv
      iswings i | i `elem` isInt = atLeastOne' $snd $ swings!!i
                | otherwise = snd $ swings!!i
      funits d i | d/=i = 0
                 | i `elem` isInt = atLeastOne' $ snd $ swings!!i
                 | otherwise = snd $ swings!!i 
      units = traceit "\n\nunits" $ flip map [0..n-1] $ \d -> L.buildVector n $ funits d
      --http://www.caspur.it/risorse/softappl/doc/sas_docs/ormp/chap5/sect28.htm
      fhess (i,j) | i>=j = 
                      ((f $ xv + units!!i + units!!j)
                       - (f $ xv + units!!i - units!!j)
                       - (f $ xv - units!!i + units!!j)
                       + (f $ xv - units!!i - units!!j) ) 
                      / (4*(iswings i) * (iswings j))
                  | otherwise = 0.0    
      hess1= L.buildMatrix n n fhess 
      hess2 = L.buildMatrix n n $ \(i,j) ->if i>=j then hess1 L.@@>(i,j) 
                                                 else hess1 L.@@>(j,i) 
  in (L.fromList (map (fst) swings), hess2)


setM i j v m = L.buildMatrix (L.rows m) (L.cols m)$ \(i',j') ->if i'==i && j==j' then v
                                                               else m L.@@>(i',j') 



atLeastOne' :: Double -> Double
atLeastOne' x | isNaN x || isInfinite x = 1.0
              | x < -1.0 || x > 1.0 = realToFrac $ round x
              | x < 0 = -1.0
              | otherwise = 1.0



--traceit s x = trace (s++": "++show x) x

setN x v = L.buildVector 6 $ \ix -> if ix == 0 then x else v L.@> ix

cut = 500

sigSlope :: Signal Double -> Double
sigSlope (Signal dt _ sv) = fst $ runStat regressF $ zip [0,dt..] $ L.toList sv

weigh _ [] = []
weigh t_hat (pt@(t_data,y):rest) 
  | abs(t_data-t_hat) > cut = weigh t_hat rest
  | otherwise = replicate (5 - (floor $ abs(t_data-t_hat)/cut * 5)) pt ++ weigh t_hat rest


weighRegression :: [(Double,Double)] -> Double -> Double
weighRegression pts x = y where
   myPts = weigh x pts -- filter (\(t,y)-> abs(x-t) < cut) pts
   (slope,offset) = runStat regressF myPts
   y = slope*x + offset

zscore :: [(Double,Double)] -> (Double,Double) -> Double
zscore tsamps (t,amp) = abs z where
  relevant = filter ((< 200) . abs . (`subtract` t) . fst) tsamps
  (slope, offset) = runStat regressF $ relevant
  detrended = map (\(t,amp)-> amp - (slope*t + offset)) relevant
  (mean, sd) = runStat meanSDF $ detrended
  z = (amp - (slope*t + offset))/ sd

localVar :: [(Double,Double)] -> (Double,Double) -> Double
localVar tsamps (t,amp) = sd*sd where
  relevant = filter ((< 200) . abs . (`subtract` t) . fst) tsamps
  (slope, offset) = runStat regressF $ relevant
  detrended = map (\(t,amp)-> amp - (slope*t + offset)) relevant
  (mean, sd) = runStat meanSDF $ detrended
  z = (amp - (slope*t + offset))/ sd
 

showNPQV' am = showNPQV (ampPar am)++" lh: "++show (lastLike am) ++ " rate " ++ show rate
  where rate = realToFrac (count_accept am) / realToFrac (count am)
showNPQV :: L.Vector Double -> String
showNPQV (L.toList -> [n, logcv, phi, logq]) 
  = intercalate "\t" $ map (minWidth 8 . accushow) [n, exp logcv, phi, exp logq]
  

minWidth n s | len < n = s ++ replicate (n - len) ' '
             | otherwise = s
  where len = length s

simq = 4.000e-4
--simn = 50

--simntrials = 1000
simt0 = 0.035


logthetaHat = 0.15 --4.630e-3
sigmaHat = 1.6 --6.260
logobsHat = 2.7e-4

hatPars = [1.597389498633494,1.142442855327917,-8.1510424681037]

fakesam simn ntrials = return 0.150>>=(\cv-> 
        (return 0.800)>>=(\phi-> 
        (return 0.200)>>=(\plo-> 
        (return (simq*(100/realToFrac simn)))>>=(\q-> 
        (return 170.000)>>=(\tc-> 
        (for 1 ntrials)$(\i-> 
            let p = phi-((phi-plo)/(1.000+(exp (soffset*sslope-sslope*realToFrac i)))) 
            in ((((binGauss simn p) q) cv) 0.000)>>=(\amp-> 
               (return (0-60.000))>>=(\vstart-> 
               ((gpByChol dt) 
                          (\t-> vstart+(amp*(((((step (t-simt0))*tc)*tc)*(t-simt0))*(exp ((0.000-(t-simt0))*tc)))))) 
                          cholm))))))))
  where cholm = chol $ mkCovM (np+1) hatPars
        soffset = ntrials / 2
        sslope = 8.000e-3 * (1000/ntrials) 

quantalEvent t = simq*(((((step (t-simt0))*170)*170)*(t-simt0))*(exp ((0.000-(t-simt0))*170)))
stagger (i, Signal dt t0 sig) = Signal dt i sig

pad (c:[]) = '0':c:[]
pad cs = cs

--fst3 (x,_,_) = x

whenContinues sess mma = do
      conts <- durations "continues" "foo"
      sessid <- getSessionName
      case conts of
        [] -> if (sess `isPrefixOf` sessid || sessid `isPrefixOf` sess) then mma else return Nothing
        (_,s):_ | s `isPrefixOf` sess -> mma
                | otherwise -> return Nothing
  


getBurnIn sess = case find (\(s,v) -> s `isPrefixOf` sess) burnIn of
                  Just (s,v) -> v
                  Nothing -> 0

getStartN sess = case find (\(s,v) -> s `isPrefixOf` sess) startNs of
                  Just (s,v) -> v
                  Nothing -> 30

getFilter sess = case find (\(s,v) -> s `isPrefixOf` sess) filters of
                  Just (s,v) -> v
                  Nothing -> const True

{-posteriorTop pcurve amparloops v = -- ((n,cv,slope,offset,phi,plo,q,tc,t0), loopvals) = 
 oneToLogPdf (800) n
 +uniformLogPdf 0.00001 0.5 cv
 +uniformLogPdf (0) (1) phi
 +uniformLogPdf (0.000) (1) q
 +(sum $ (flip map) (zip pcurve  amps) $ \(pcurveVal, amp)->
    let p=pcurveVal * phi 
        --nr::Int = round $ relfrac *p* realToFrac n
    in binGaussFrom1LogPdf n p q cv 0 amp) --  + -- binomialLogProb (n) (p) nr +
       --normalLogPdf (realToFrac nr*q) (varToTau (q*cv*q*cv*(realToFrac nr))) amp )
  where n = round $ v @> 0
        cv = exp $ v @> 1
        phi = v @> 2
        q = exp $ v @> 3
        --relfracs = map ( (@>0) .  ampPar) amparloops
        amps =map ( (@>0) .  ampPar) amparloops

-}
getSess def = do
  args <- getArgs 
  case args of 
    [] -> return def
    s : _ -> return s

datasess = 
 --words "00c9bd 0ca3a9 84b41c 22b152 512f48 7b8f60 b34863 b62b8f cf96ab fcb952 57246a"
 words "00c9bd 0ca3a9 84b41c 57246a 22b152 512f48 7b8f60 fcb952 b34863 cf96ab b62b8f"  -- fcb952 b62b8f

burnIn = [("22b", 8000), ("b62", 6000), ("cf96", 4000), 
          ("b34", 3000), ("00c9",10000), ("84", 15000), ("572", 13000)]

slopeFilter = [("0ca", 4)]

startNs = [("00c9", 50), ("84", 40), ("512", 50) ]

filters = [("b34", \(t,amp)-> t>5300), 
           ("cf96", \(t,amp)-> t>400),
           ("84b", (>500) . fst)]

nsol = 4

getWf sess = do
  nms <- fmap read $  readFile (take 6 sess++"/sessions")
  sigs' <- fmap concat $ forM nms $ \sessNm-> do 
            LoadSignals sigs <- decodeFile $ take 6 sess++"/sigs_"++take 6 sessNm++"_epsps" 
            return sigs
  let sigs = case lookup (take 3 sess) slopeFilter of
                Nothing -> sigs'
                Just x -> filter ((<x) . abs . sigSlope) sigs'
  let wf@(Signal _ _ sv) = zeroSigInit 0.03 $ baselineSigFrom 0.02 0.032 $ averageSigs $ sigs
  return $ (wf, foldl1' max $ L.toList sv, sigs')

getWf' sess = do
  nms <- fmap read $  readFile ("sessions")
  sigs' <- fmap concat $ forM nms $ \sessNm-> do 
            LoadSignals sigs <- decodeFile $ "sigs_"++take 6 sessNm++"_epsps" 
            return sigs
  let sigs = case lookup (take 3 sess) slopeFilter of
                Nothing -> sigs'
                Just x -> filter ((<x) . abs . sigSlope) sigs'
  let wf@(Signal _ _ sv) = zeroSigInit 0.03 $ baselineSigFrom 0.02  0.032 $ averageSigs $ sigs
  return $ (wf, foldl1' max $ L.toList sv, sigs')




{-posteriorLoop wf invDetails ampartop pcurveVal sig v 
  = binomialLogProb (n) (p) nr +
    normalLogPdf (realToFrac nr*q) (varToTau (q*cv*q*cv*(realToFrac nr))) amp +
    ouSynapseLogPdf invDetails (scaleSig (v0) amp wf) sig
 where v0 = v@> 2
       amp = v@>1       
       nr = round $ v@>0
       vt = ampPar ampartop
       n = round $ vt @> 0
       cv = exp $ vt @> 1
       phi = vt @> 2
       q = exp $ vt @> 3
       p = pcurveVal * phi
 
posteriorLoop' sd amtop pcurveVal sigAmpMean v 
  = binGaussFrom1LogPdf n p q cv 0 amp + --binomialLogProb' (n) (p) nr +
--    normalLogPdf (realToFrac nr*q) (varToTau (q*cv*q*cv*(realToFrac nr))) amp +
    --ouSynapseLogPdf invDetails (scaleSig (v0) amp wf) sig
    normalLogPdf sigAmpMean (varToTau $ sd*sd ) amp
 where --v0 = v@> 2
       amp = v@>0       
       --relfrac =  v@>0
       --nr = round $ relfrac *p* realToFrac n
       vtop = ampPar amtop
--       vQCV = ampPar amtopQCV
       n = round $ vtop @> 0
       cv = exp $ vtop @> 1
       phi = vtop @> 2
       q = exp $ vtop @> 3
       p = pcurveVal * phi
       
updateG wf invDetails pcurve sigs (ampartop,amparloops) = do
   newtop <- adaMetNoCacheP False (posteriorTop pcurve amparloops) ampartop
   newloops <- forM (zip3 pcurve sigs amparloops) $ \(pcurveVal, sig, ampar) ->
                 adaMetNoCacheP False 
                        (posteriorLoop wf invDetails newtop pcurveVal sig) 
                        ampar
   return (newtop, newloops)

updateG' sd amps pcurve topcool (amtop, amparloops) = do
   newtop <- adaMetNoCacheP False ((/topcool) . posteriorTop pcurve amparloops) amtop
   --newtopQCV <- adaMetNoCacheP False ((/topcool) . posteriorTopQCV pcurve newtopNP amparloops) amtopQCV
   newloops <- forM (zip3 pcurve amps amparloops) $ \(pcurveVal, amp, ampar) ->
                 adaMetNoCacheP False 
                        (posteriorLoop' sd newtop pcurveVal amp ) 
                        ampar
   return (newtop, newloops)


runGibbs 0 wf invDetails pcurve sigs (ampartop,amparloops) xs = return $ reverse xs
runGibbs n wf invDetails pcurve sigs (ampartop,amparloops) xs = do
      (newtop, newloops)<- updateG wf invDetails pcurve sigs (ampartop,amparloops)
      runGibbs (n-1) wf invDetails pcurve sigs (newtop, newloops) (ampPar newtop:xs)

runGibbs' [] sd amps pcurve pars xs = return $ reverse xs
runGibbs' ((0,_):rest) sd amps pcurve pars xs = do
      runGibbs' (rest) sd amps pcurve pars xs
runGibbs' ((n,topcool):rest) sd amps pcurve pars xs = do
      newpars@(p1,p3)<- updateG' sd amps pcurve topcool pars
      runGibbs' ((n-2,topcool):rest) sd amps pcurve newpars
                (L.join (map ampPar [p1,head p3]):xs)

-}
reset_counts n ampar = ampar { count = n, count_accept = n `div` 2} 

shrink x ampar = ampar {ampCov = L.scale (recip x) $ ampCov ampar} 

wToRGB w 
 | w >= 380 && w < 440 = 
    (-(w - 440) / (440 - 380)
    ,0.0
    ,1.0)
 | w >= 440 && w < 490 = 
    (0.0
    ,(w - 440) / (490 - 440)
    ,1.0)
 | w >= 490 && w < 510 = 
    (0.0
    ,1.0
    ,-(w - 510) / (510 - 490))
 | w >= 510 && w < 580 = 
    ((w - 510) / (580 - 510)
    ,1.0
    ,0.0)
 | w >= 580 && w < 645 = 
    (1.0
    ,-(w - 645) / (645 - 580)
    ,0.0)
 | w >= 645 && w <= 780 = 
    (1.0
    ,0.0
    ,0.0)
 | otherwise = 
    (0.0
    ,0.0
    ,0.0)

{-posteriorTopNP pcurve amTopQCV amparloops v = -- ((n,cv,slope,offset,phi,plo,q,tc,t0), loopvals) = 
 oneToLogPdf (800) n
 +uniformLogPdf (0) (1) phi
 +(sum $ (flip map) (zip3 pcurve relfracs amps) $ \(pcurveVal, relfrac, amp)->
    let p=pcurveVal * phi 
        nr::Int = round $ relfrac *p* realToFrac n
    in binomialLogProb' (n) (p) nr)
  where n = round $ v @> 0
        --cv = exp $ v @> 1
        phi = v @> 1
        --q = exp $ v @> 3
        relfracs = map ( (@>0) .  ampPar) amparloops
        amps =map ( (@>1) .  ampPar) amparloops

posteriorTopQCV pcurve amTopNP amparloops v = -- ((n,cv,slope,offset,phi,plo,q,tc,t0), loopvals) = 
 normalLogPdf (0.2) (varToTau $ 0.05*0.05)cv
 +uniformLogPdf (0.000) (1) q
 +sm
  where n = round $ (ampPar amTopNP) @> 0
        cv = exp $ v @> 0
        phi = (ampPar amTopNP) @> 1
        q = exp $ v @> 1
        relfracs = map ( (@>0) .  ampPar) amparloops
        amps =map ( (@>1) .  ampPar) amparloops
        sm = (sum $ (flip map) (zip3 pcurve relfracs amps) $ \(pcurveVal, relfrac, amp)->
                let p=pcurveVal * phi 
                    nr::Int = round $ relfrac *p* realToFrac n
                    fakecv= 0.2
                    it = normalLogPdf (realToFrac nr*q) (varToTau (q*fakecv*q*fakecv*(realToFrac nr))) amp 
                in it ) -- `seq` if isNaN it then trace (show (p, relfrac, amp, nr)) it else it)
-}