function INP(N)

Pos = 100.0 * lhsdesign(N,2);
Den = rand(N,1);

LHS = [Pos, Den];

str = sprintf('Particle%d.dat',N);
dlmwrite(str,LHS,'delimiter',' ','precision',5)