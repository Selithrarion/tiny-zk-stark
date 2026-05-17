### tiny-zk-stark

trace -> AIR -> FRI -> proof  

a study, building the core components of STARK from scratch (then replaced some stuff with arkworks)
* ntt/intt (fast fourier transform but with integers)
* air (frontend)
* fri (backend)
* poseidon, babybear, merkle trees
* polynomial logic

todo:  
* lookup tables
* optimize poly module or replace with arkworks
* use rayon?

poseidon2 with babybear params and merkle tree taken from here https://github.com/HorizenLabs/poseidon2/tree/main