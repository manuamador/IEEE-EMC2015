const c=299792458.
const mu0=4*pi*1e-7
const mur=1
const eps0=1./(mu0*c^2)
const epsr=1

using PyPlot

function urc(Q,f,f0)
    w=f/f0
    I=1/sqrt(1+Q^2*(w-1/w)^2)
    return I
end

const a=2.9#3.11
const b=3.7#3.62
const d=8.7#4.84
const f0=100e6
const f1=1200e6

#rf=.01
#nf=int(1/log(1+rf)*log(f1/f0))
#f=zeros(nf)
#for i=1:nf
#  f[i]=f0*(1+rf)^i
#end


const df=.1e6
const nf=int((f1-f0)/df)+1
const f=linspace(f0,f1,nf)#array([20e6,40e6,45e6,56e6,65e6])#linspace(20e6,.3e9,281)#linspace(20e6,3e8,1121)
const V=a*b*d
const Q4=16*pi^2*V*f.^3/(1*c^3)
const Q1=3*V./(2*1*(2./sqrt(2*pi*f.*mur*mu0*1*1e5)))/200
const Q=1./(1./Q1+1./Q4)

modes=zeros(nf)


const Npt=50

const x=rand(Npt)*(a-1)+0.5
const y=rand(Npt)*(b-1)+0.5
const z=rand(Npt)*(d-1)+0.5

E0=1
H0=1/(c*mu0*mur)
Ex=zeros(Complex128,Npt,nf)
Ey=zeros(Complex128,Npt,nf)
Ez=zeros(Complex128,Npt,nf)
for u=1:nf
  n_x=int(f[u]*(1+100/Q[u])*2*a/c)
  n_y=int(f[u]*(1+100/Q[u])*2*b/c)
  n_z=int(f[u]*(1+100/Q[u])*2*d/c)
  M=[]
  N=[]
  P=[]
  F=[]
  AMP=[]
  for mm=0:n_x
    for nn=0:n_y
      for pp=0:n_z
        fmnp=c/2*sqrt((mm/a)^2+(nn/b)^2+(pp/d)^2)
        if abs(fmnp-f[u])<100*f[u]/Q[u]
          Amp=urc(Q[u],f[u],fmnp)
          M=vcat(M,mm)
          N=vcat(N,nn)
          P=vcat(P,pp)
          F=vcat(F,fmnp)
          AMP=vcat(AMP,Amp)
        end
      end
    end
  end

  nm=length(M)
  println("$(round(f[u]*100/1e6)/100) MHz, $nm resonant freq.")
  for k=1:nm
    fmnp=F[k]
    Amp=AMP[k]
    kx=M[k]*pi/a
    ky=N[k]*pi/b
    kz=P[k]*pi/d
    k2mnp=kx^2+ky^2+kz^2
    for i=1:Npt
      ETMx_x=0
      ETMx_y=0
      ETMx_z=0
      ETMy_x=0
      ETMy_y=0
      ETMy_z=0
      ETMz_x=0
      ETMz_y=0
      ETMz_z=0
      ETEx_y=0
      ETEx_z=0
      ETEy_x=0
      ETEy_z=0
      ETEz_y=0
      ETEz_x=0
      ssc_z=sin(kx*x[i])*sin(ky*y[i])*cos(kz*z[i])
      scs_z=sin(kx*x[i])*cos(ky*y[i])*sin(kz*z[i])
      css_z=cos(kx*x[i])*sin(ky*y[i])*sin(kz*z[i])

      ssc_x=sin(ky*y[i])*sin(kz*z[i])*cos(kx*x[i])
      scs_x=sin(ky*y[i])*cos(kz*z[i])*sin(kx*x[i])
      css_x=cos(ky*y[i])*sin(kz*z[i])*sin(kx*x[i])

      ssc_y=sin(kz*z[i])*sin(kx*x[i])*cos(ky*y[i])
      scs_y=sin(kz*z[i])*cos(kx*x[i])*sin(ky*y[i])
      css_y=cos(kz*z[i])*sin(kx*x[i])*sin(ky*y[i])

      if M[k]>0 && N[k]>0 && P[k]>0
        #TM modes z
        ETMx_z = -E0*Amp*kx*kz/(k2mnp-kz^2)*css_z
        ETMy_z = E0*Amp*ky*kz/(k2mnp-kz^2)*scs_z
        ETMz_z = E0*Amp*ssc_z
        #TE modes z
        ETEx_z = -1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kz^2)*css_z
        ETEy_z = 1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-kz^2)*scs_z

        #TM modes x
        ETMy_x = -E0*Amp*ky*kx/(k2mnp-kx^2)*css_x
        ETMz_x = E0*Amp*kz*kx/(k2mnp-kx^2)*scs_x
        ETMx_x = E0*Amp*ssc_x
        #TE modes x
        ETEy_x = -1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-kx^2)*css_x
        ETEz_x = 1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kx^2)*scs_x

        ##TM modes y
        ETMz_y = -E0*Amp*kz*ky/(k2mnp-ky^2)*css_y
        ETMx_y = E0*Amp*kx*ky/(k2mnp-ky^2)*scs_y
        ETMy_y = E0*Amp*ssc_y
        ##TE modes y
        ETEz_y = -1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-ky^2)*css_y
        ETEx_y = 1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-ky^2)*scs_y
      end

      if M[k]>0 && N[k]>0 && P[k]==0
        ETMz_z = E0*Amp*ssc_z #TM mode z
        ETEz_y = -1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-ky^2)*css_y #TE mode y
        ETEz_x = 1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kx^2)*scs_x #TE mode x
      end
      if M[k]==0 && N[k]>0 && P[k]>0
        ETEx_z = -1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kz^2)*css_z #TE mode z
        ETMx_x = E0*Amp*ssc_x #TM mode x
        ETEx_y = 1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-ky^2)*scs_y #TE mode y
      end
      if M[k]>0 && N[k]==0 && P[k]>0
        ETEy_z = 1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-kz^2)*scs_z #TE mode z
        ETMy_y = E0*Amp*ssc_y #TM mode y
        ETEy_x = -1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-kx^2)*css_x #TE mode x
      end
      Ex[i,u]+=ETMx_x+ETMx_y+ETMx_z+ETEx_y+ETEx_z
      Ey[i,u]+=ETMy_x+ETMy_y+ETMy_z+ETEy_x+ETEy_z
      Ez[i,u]+=ETMz_x+ETMz_y+ETMz_z+ETEz_y+ETEz_x
    end
  end
end

using NPZ
npzwrite("E_fullmodes50.npz", ["Ex" => Ex,"Ey" => Ey,"Ez" => Ez, "f" => f, "Q" => Q, "x" => x,"y" => y,"z" => z])

#npzwrite("E_fullmodes50.npz", ["Ex" => Ex, "f" => f, "Q" => Q, "x" => x,"y" => y,"z" => z])

function CDF(X)
  x=sort(X)
  y=cumsum(x)/sum(x)
  return y
end


function ADTest(X,distribution,alpha=0.05)
"""
Anderson Darling statistical test for Rayleigh or exponential distribution
return 0 if the distribution follows the given distribution, 1 otherwise.

X: vector, N samples (30<N<150)
distribution: string 'Rayleigh' or 'Exponential'
alpha: typical statistical error
"""

	N=length(X)

	#Stephen's critical values
	if alpha==0.15
		CV=0.922
	elseif alpha==0.10
		CV=1.078
	elseif alpha==0.05
		CV=1.341
	elseif alpha==0.025
		CV=1.606
	elseif alpha==0.01
		CV=1.957
	end
	if distribution=="Rayleigh"
		b=mean(X)/sqrt(pi/2); #get the shape parameter
		x= minimum(X):(maximum(X)-minimum(X))/(N-1):maximum(X);
		y=sort(X);
		CDF = 1.-exp(-y.^2/2b^2); # theoretical CDF
		Stat=log(CDF)+log(1.-flipud(CDF));
	elseif distribution=="Exponential"
		b=1/mean(X); #get the scale parameter
		x= minimum(X):(maximum(X)-minimum(X))/(N-1):maximum(X);
		y=sort(X);
		CDF = 1.-exp(-b.*y); # theoretical CDF
		Stat=log(CDF)+log(1.-flipud(CDF));
	end

	A=[];
    for i=1:N
        A=[A;(2*i-1)*Stat[i]];
    end
    Astat=-sum(A)/N-N;
    if Astat*(1+0.6/N)>CV #test
        result=1;
    else
        result=0;
    end
  return result,Astat
end

ADresultx=zeros(nf)
ADresulty=zeros(nf)
ADresultz=zeros(nf)


for u=1:nf
  ADresultx[u]=ADTest(abs(Ex[:,u])/maximum(abs(Ex[:,u])),"Rayleigh",0.05)[1]
  ADresulty[u]=ADTest(abs(Ey[:,u])/maximum(abs(Ey[:,u])),"Rayleigh",0.05)[1]
  ADresultz[u]=ADTest(abs(Ez[:,u])/maximum(abs(Ez[:,u])),"Rayleigh",0.05)[1]
end

#moving average
w=100
avADx=zeros(nf-w)*NaN
avADy=zeros(nf-w)*NaN
avADz=zeros(nf-w)*NaN


for u=1:nf-w
  avADx[u]=mean(ADresultx[u:u+w])
  avADy[u]=mean(ADresulty[u:u+w])
  avADz[u]=mean(ADresultz[u:u+w])
end

fch=linspace(200,1100,11)
rrch=[0.82,0.58,0.70,0.45,0.29,0.1,0.2,0.15,0.15,0.13,0.07]


figure(3)
plot(f[1:end-w]/1e6,avADx,".b",label="\$E_x\$")
plot(f[1:end-w]/1e6,avADy,"+g",label="\$E_y\$")
plot(f[1:end-w]/1e6,avADz,"xr",label="\$E_z\$")
plot(fch,rrch,color="k","s",markerfacecolor="w", markersize=10,label="\$E_R\$, Meas.")
grid()
xlabel("\$ f\$/MHz")
ylabel("Rejection rate")
legend()
savefig("fig3.pdf",bbox_inches="tight")
