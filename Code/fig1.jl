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
const f1=2000e6
rf=.005
nf=int(1/log(1+rf)*log(f1/f0))
f=zeros(nf)
for i=1:nf
  f[i]=f0*(1+rf)^i
end

const V=a*b*d
const Q4=16*pi^2*V*f.^3/(1*c^3)
const Q1=3*V./(2*1*(2./sqrt(2*pi*f.*mur*mu0*1*1e5)))/200
const Q=1./(1./Q1+1./Q4)
modes=zeros(nf)
const x=1
const y=2
const z=1

E0=1
H0=1/(c*mu0*mur)
Ex=zeros(Complex128,nf)
Ey=zeros(Complex128,nf)
Ez=zeros(Complex128,nf)
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
        Amp=urc(Q[u],f[u],fmnp)
        if Amp>0.005
          M=vcat(M,mm)
          N=vcat(N,nn)
          P=vcat(P,pp)
          F=vcat(F,fmnp)
          AMP=vcat(AMP,Amp)
        end
      end
    end
  end
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
  nm=length(M)
  println("$(round(f[u]*100/1e6)/100) MHz, $nm resonant freq.")
  for k=1:nm
    fmnp=F[k]
    Amp=AMP[k]
    #kx=M[k]*pi/a
    #ky=N[k]*pi/b
    #kz=P[k]*pi/d
    #k2mnp=kx^2+ky^2+kz^2
    #ssc_z=sin(kx*x)*sin(ky*y)*cos(kz*z)
    #scs_z=sin(kx*x)*cos(ky*y)*sin(kz*z)
    #css_z=cos(kx*x)*sin(ky*y)*sin(kz*z)

    #ssc_x=sin(ky*y)*sin(kz*z)*cos(kx*x)
    #scs_x=sin(ky*y)*cos(kz*z)*sin(kx*x)
    #css_x=cos(ky*y)*sin(kz*z)*sin(kx*x)

    #ssc_y=sin(kz*z)*sin(kx*x)*cos(ky*y)
    #scs_y=sin(kz*z)*cos(kx*x)*sin(ky*y)
    #css_y=cos(kz*z)*sin(kx*x)*sin(ky*y)

    if M[k]>0 && N[k]>0 && P[k]>0
      #TM modes z
      #ETMx_z += -E0*Amp*kx*kz/(k2mnp-kz^2)*css_z
      #ETMy_z += E0*Amp*ky*kz/(k2mnp-kz^2)*scs_z
      #ETMz_z += E0*Amp*ssc_z
      ##TE modes z
      #ETEx_z += -1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kz^2)*css_z
      #ETEy_z += 1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-kz^2)*scs_z
      modes[u] += 2

      ##TM modes x
      #ETMy_x += -E0*Amp*ky*kx/(k2mnp-kx^2)*css_x
      #ETMz_x += E0*Amp*kz*kx/(k2mnp-kx^2)*scs_x
      #ETMx_x += E0*Amp*ssc_x
      #TE modes x
      #ETEy_x += -1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-kx^2)*css_x
      #ETEz_x += 1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kx^2)*scs_x
      modes[u] += 2

      ##TM modes y
      #ETMz_y += -E0*Amp*kz*ky/(k2mnp-ky^2)*css_y
      #ETMx_y += E0*Amp*kx*ky/(k2mnp-ky^2)*scs_y
      #ETMy_y += E0*Amp*ssc_y
      ##TE modes y
      #ETEz_y += -1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-ky^2)*css_y
      #ETEx_y += 1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-ky^2)*scs_y
      modes[u] += 2
    end

    if M[k]>0 && N[k]>0 && P[k]==0
      #ETMz_z+=E0*Amp*ssc_z #TM mode z
      #ETEz_y+=-1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-ky^2)*css_y #TE mode y
      #ETEz_x+=1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kx^2)*scs_x #TE mode x
      modes[u]+=3
    end
    if M[k]==0 && N[k]>0 && P[k]>0
      #ETEx_z+=-1im*2*pi*fmnp*mur*mu0*ky*H0*Amp/(k2mnp-kz^2)*css_z #TE mode z
      #ETMx_x+=E0*Amp*ssc_x #TM mode x
      #ETEx_y+=1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-ky^2)*scs_y #TE mode y
      modes[u]+=3
    end
    if M[k]>0 && N[k]==0 && P[k]>0
      #ETEy_z+=1im*2*pi*fmnp*mur*mu0*kx*H0*Amp/(k2mnp-kz^2)*scs_z #TE mode z
      #ETMy_y+=E0*Amp*ssc_y #TM mode y
      #ETEy_x+=-1im*2*pi*fmnp*mur*mu0*kz*H0*Amp/(k2mnp-kx^2)*css_x #TE mode x
      modes[u]+=3
    end
  end
  #Ex[u]=ETMx_x+ETMx_y+ETMx_z+ETEx_y+ETEx_z
  #Ey[u]=ETMy_x+ETMy_y+ETMy_z+ETEy_x+ETEy_z
  #Ez[u]=ETMz_x+ETMz_y+ETMz_z+ETEz_y+ETEz_x
end

figure(num=1)#,figsize=(6,4.5))
subplot(211)
title("\$ \\Delta f \$")
loglog(f,100*f./Q)#,label="\$ \\vert E_x\\vert\$")
ylabel("Hz")
xlim(f0,f1)
grid(which="both")
subplot(212)
title("Number of modes")
loglog(f,modes)#,label="\$ \\vert E_x\\vert\$")
#ylabel("A.U")
xlim(f0,f1)
ylim(10,100000)
grid(which="both")
xlabel("\$f\$ in Hz")
savefig("fig1.pdf",bbox_inches="tight")

