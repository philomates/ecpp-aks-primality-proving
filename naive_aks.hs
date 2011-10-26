-- Phillip Mates
-- Attempt at naive AKS
-- Work in Progress

aks :: Integer -> Bool
aks n
  | powerOfPrime n 2 1 == True = False
  | otherwise =
    case findR 2 n of
      Just r -> fermatsCheck 1 r n
      Nothing -> False

findR :: Integer -> Integer -> Maybe Integer
findR r n
  | r >= n = Just r
  | gcd n r /= 1 = Nothing
  | isPrime r =
    let q = 3 -- FIXME: largest factor of r-1
    in
      if and [q > 4 * floor (sqrt (fromIntegral r) * logBase 2 (fromIntegral n))
            , n^(r-1 `div` q) `mod` r /= 1]
        then Just r
        else findR (r+1) n

  | otherwise = findR (r+1) n

-- b^k (mod m)
-- Via paper titled "Implementing Public-Key Cryptography in Haskell"
expm :: Integer -> Integer -> Integer -> Integer
expm m b k =
  let
    ex a k s
      | k == 0 = s
      | k `mod` 2 == 0 = ((ex (a*a `mod` m)) (k `div` 2)) s
      | otherwise = ((ex (a*a `mod` m)) (k `div` 2)) (s*a `mod` m)
  in ex b k 1

-- naive slow b^k (mod m)
expm' :: Integer -> Integer -> Integer -> Integer
expm' m b k = (b ^ k) `mod` m


-- TODO: primality check for smaller numbers
isPrime :: Integer -> Bool
isPrime = undefined

fermatsCheck :: Integer -> Integer -> Integer -> Bool
fermatsCheck a r n
  | a < 2 * (floorSqrt r) * floorLn n =
    if a^n `mod` n == a -- FIXME: incorrect check!
       then False
       else fermatsCheck (a+1) r n
  | otherwise = True

powerOfPrime :: Integer -> Integer -> Integer -> Bool
powerOfPrime n a b
  | and [a < floorSqrt n
       , b > flooLn n] = powerOfPrime n (a+1) 1
  | and [a >= floorSqrt n
       , b > flooLn n] = False
  | otherwise =
    case compare n (a^b) of
      EQ -> True
      LT -> powerOfPrime n a (b+1)
      GT -> powerOfPrime n (a+1) 2

floorSqrt :: Integer -> Integer
floorSqrt n = floor $ sqrt $ fromIntegral n

floorLn :: Integer -> Integer
floorLn n = floor $ logBase 2 (fromIntegral n)
