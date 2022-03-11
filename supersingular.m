num_prime_factors := 80;
num_primes_used_in_walk := 12;
length_walk := 2^5;
num_tests := 20;

for t in [1..num_tests] do

  // generate test curve: csidh curve over Fp with p = 1 mod 4
  // and p = 2 mod 3 (automatic if num_prime_factors > 1)

  primes := [];
  ell := 1; prod := 1;
  for i in [1..num_prime_factors-1] do
    ell := NextPrime(ell);
    primes cat:= [ell];
    prod *:= ell;
  end for;
  ell := Ceiling(Random([10..100])*ell/Random([1..10]));
  repeat
    ell := NextPrime(ell);
  until IsPrime(prod*ell - 1);
  primes cat:= [ell];
  p := prod*ell - 1; // note: automatically p = 1 mod 4 and p = 2 mod 3
  Fp := GF(p);
  R<x> := PolynomialRing(Fp);
  E := EllipticCurve(x^3 + 1);

  // take random CSIDH walk and keep track of character value

  Enew := E;
  chi := 1;
  walkdone := [];
  for i in [1..length_walk] do
    ell := Random(primes[2..num_primes_used_in_walk + 1]); // skip ell = 2
    cofac := (p + 1) div ell;
    if Random(1) eq 1 then
      twist := true;
      Enew := QuadraticTwist(Enew); // note: no issues with twist if p = 1 mod 4
    else
      twist := false;
    end if;
    repeat
      P := cofac*Random(Enew);
    until P ne Enew ! 0;
    kerpol := &*[x - (i*P)[1] : i in [1..ell div 2]];
    Enew := IsogenyFromKernel(Enew, kerpol);
    if twist then
      Enew := QuadraticTwist(Enew);
    end if;
    chi *:= (-1)^((ell - 1) div 2);
  end for;

  // compute Weil pairings on E and Enew, use sigma = frob

  s := 2; r := 2*s;
  Fq := GF(p, r); // 4-torsion is defined over Fq
  onebranch := p^s - (-1)^s^2; // square root of number of pts
  cofac := onebranch div 4;

  E := BaseChange(E, Fq);
  Enew := BaseChange(Enew, Fq);

  function frob(P, p)
    return Curve(P) ! [P[1]^p, P[2]^p, P[3]^p];
  end function;

  repeat
    P := cofac*Random(E);
    Pfrob := frob(P, p);
  until 2*P ne 2*Pfrob;
  mu := WeilPairing(P, Pfrob, 4);

  repeat
    P := cofac*Random(Enew);
    Pfrob := frob(P, p);
  until 2*P ne 2*Pfrob;
  munew := WeilPairing(P, Pfrob, 4);

  print "Bit length of field size:", Round(Log(2, p));
  print "Is mu_E primitive?", mu^2 ne 1;
  print "Is mu_Enew primitive?", munew^2 ne 1;
  test := munew eq mu^chi;
  print "Is mu_Enew = mu_E^chi?", test, "( chi was", chi, ")";
  print "-----";

  if not test then // if we trust result, then test not needed :-)
    print "Problem with", E, Enew;
    break;
  end if;
end for;
