import sys
print(sys.argv[1]);
head='''@RULE 1dca
@TABLE
n_states:3
neighborhood:Moore
symmetries:none
var m=0
var p=2
var a={m,1}
var b={m,1}
var c={m,1}
var d={m,1}
var e={m,1}
var f={m,1}
var g={m,1}
var h={m,1}

m,p,p,m,m,m,f,g,p,%s
m,p,p,m,m,m,f,g,1,%s
m,p,1,m,m,m,f,g,p,%s
m,p,1,m,m,m,f,g,1,%s
m,1,p,m,m,m,f,g,p,%s
m,1,p,m,m,m,f,g,1,%s
m,1,1,m,m,m,f,g,p,%s
m,1,1,m,m,m,f,g,1,%s

1,a,b,d,e,f,g,h,c,1
p,a,b,d,e,f,g,h,c,p

''';
file='1dca.rule';
rnum=int(sys.argv[1]);
r=bin(rnum);
r=r[:1:-1];
print(r)
r+='0'*(8-len(r));
r=r[::-1];
rule=r;
# rule=[0, 1, 1, 1, 1, 1, 0, 0];
# dct={0:'p',1:'1'};
dct=['p','1'];
with open(file,'w') as f:
	f.write(head%tuple([dct[int(x)] for x in rule]));
	# for i in rule:
		# f.write(dct[i]);