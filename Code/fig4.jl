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

const df=2e6
const nf=1
const f=[274e6]#[224e6]#[208e6]#[382e6]#linspace(f0,f1,nf)#array([20e6,40e6,45e6,56e6,65e6])#linspace(20e6,.3e9,281)#linspace(20e6,3e8,1121)
const V=a*b*d
const Q4=16*pi^2*V*f.^3/(1*c^3)
const Q1=3*V./(2*1*(2./sqrt(2*pi*f.*mur*mu0*1*1e5)))/500
const Q=1./(1./Q1+1./Q4)


const x=1.0
dy=0.05
dz=0.05
ny=int(b/dy)+1
nz=int(d/dz)+1
const y=linspace(0,b,ny)
const z=linspace(0,d,nz)




z_b=[2,2.8,4.5]
x_b=[1,1,1]
y_b=[0.5,3,1]

n_steps=length(x_b)
#plot3D(x_pi[1:0.1:n_pos],y_pi[1:0.1:n_pos],z_pi[1:0.1:n_pos])

E0=1
H0=1/(c*mu0*mur)

Ex=zeros(Complex128,nf,n_steps,ny,nz)
#Ey=zeros(Complex128,nf,n_steps,Npt)
#Ez=zeros(Complex128,nf,n_steps,Npt)
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
        #if abs(fmnp-f[u])<100*f[u]/Q[u]
        Amp=urc(Q[u],f[u],fmnp)
        if Amp>0.01
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
  println("$nm resonant freq.")
  f_coeffTMx_x=zeros(n_steps,nm)
  f_coeffTMx_y=zeros(n_steps,nm)
  f_coeffTMx_z=zeros(n_steps,nm)
  f_coeffTEx_x=zeros(Complex128,n_steps,nm)
  f_coeffTEx_y=zeros(Complex128,n_steps,nm)
  f_coeffTEx_z=zeros(Complex128,n_steps,nm)
  for pos=1:n_steps
    ETMx_x=zeros(ny,nz)
    ETMx_y=zeros(ny,nz)
    ETMx_z=zeros(ny,nz)
    ETEx_y=zeros(Complex128,ny,nz)
    ETEx_z=zeros(Complex128,ny,nz)
    xe=x_b[pos]
    ye=y_b[pos]
    ze=z_b[pos]
    println(" Step $pos/$n_steps")
    for k=1:nm
      fmnp=F[k]
      Amp=AMP[k]
      kx=M[k]*pi/a
      ky=N[k]*pi/b
      kz=P[k]*pi/d
      k2mnp=kx^2+ky^2+kz^2
      ssc_ze=sum(sin(kx*xe).*sin(ky*ye).*cos(kz*ze))
      scs_ze=sum(sin(kx*xe).*cos(ky*ye).*sin(kz*ze))
      css_ze=sum(cos(kx*xe).*sin(ky*ye).*sin(kz*ze))
      ssc_xe=sum(sin(ky*ye).*sin(kz*ze).*cos(kx*xe))
      scs_xe=sum(sin(ky*ye).*cos(kz*ze).*sin(kx*xe))
      css_xe=sum(cos(ky*ye).*sin(kz*ze).*sin(kx*xe))
      ssc_ye=sum(sin(kz*ze).*sin(kx*xe).*cos(ky*ye))
      scs_ye=sum(sin(kz*ze).*cos(kx*xe).*sin(ky*ye))
      css_ye=sum(cos(kz*ze).*sin(kx*xe).*sin(ky*ye))

      if M[k]>0 && N[k]>0 && P[k]>0
        #TM modes z
        f_coeffTMx_z[pos,k] = -kx*kz/(k2mnp-kz^2)*css_ze
        #f_coeffTMy_z[pos,k] = ky*kz/(k2mnp-kz^2)*scs_ze
        #f_coeffTMz_z[pos,k] = ssc_ze
        #TE modes z
        f_coeffTEx_z[pos,k] = -1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze
        #f_coeffTEy_z[pos,k] = 1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze

        #TM modes x
        #f_coeffTMy_x[pos,k] = -ky*kx/(k2mnp-kx^2)*css_xe
        #f_coeffTMz_x[pos,k] = kz*kx/(k2mnp-kx^2)*scs_xe
        f_coeffTMx_x[pos,k] = ssc_xe
        #TE modes x
        #f_coeffTEy_x[pos,k] = -1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe
        #f_coeffTEz_x[pos,k] = 1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe

        ##TM modes y
        #f_coeffTMz_y[pos,k] += -kz*ky/(k2mnp-ky^2)*css_ye
        f_coeffTMx_y[pos,k] = kx*ky/(k2mnp-ky^2)*scs_ye
        #f_coeffTMy_y[pos,k] = ssc_ye
        ##TE modes y
        #f_coeffTEz_y[pos,k] = -1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye
        f_coeffTEx_y[pos,k] = 1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye
      end


      #if M[k]>0 && N[k]>0 && P[k]==0
        #f_coeffTMz_z[pos,k] = ssc_ze #TM mode z
        #f_coeffTEz_y[pos,k] = -1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye #TE mode y
        #f_coeffTEz_x[pos,k] = 1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe #TE mode x
      #end

      if M[k]==0 && N[k]>0 && P[k]>0
        f_coeffTEx_z[pos,k] += -1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze #TE mode z
        f_coeffTMx_x[pos,k] += ssc_xe #TM mode x
        f_coeffTEx_y[pos,k] += 1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye #TE mode y
      end

      #if M[k]>0 && N[k]==0 && P[k]>0
        #f_coeffTEy_z[pos,k] = 1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze #TE mode z
        #f_coeffTMy_y[pos,k] = ssc_ye #TM mode y
        #f_coeffTEy_x[pos,k] = -1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe #TE mode x
      #end

      for i=1:ny
        for j=1:nz
          ssc_z=sin(kx*x)*sin(ky*y[i])*cos(kz*z[j])
          scs_z=sin(kx*x)*cos(ky*y[i])*sin(kz*z[j])
          css_z=cos(kx*x)*sin(ky*y[i])*sin(kz*z[j])
          ssc_x=sin(ky*y[i])*sin(kz*z[j])*cos(kx*x)
          scs_x=sin(ky*y[i])*cos(kz*z[j])*sin(kx*x)
          css_x=cos(ky*y[i])*sin(kz*z[j])*sin(kx*x)
          ssc_y=sin(kz*z[j])*sin(kx*x)*cos(ky*y[i])
          scs_y=sin(kz*z[j])*cos(kx*x)*sin(ky*y[i])
          css_y=cos(kz*z[j])*sin(kx*x)*sin(ky*y[i])



          if M[k]>0 && N[k]>0 && P[k]>0#TM modes z
            #TM modes z
            ETMx_z[i,j] += E0*Amp*f_coeffTMx_z[pos,k]*css_z
            #ETMy_z[i,j] += E0*Amp*f_coeffTMy_z[pos,k]*scs_z
            #ETMz_z[i,j] += E0*Amp*f_coeffTMz_z[pos,k]*ssc_z

            #TE modes z
            ETEx_z[i,j] += H0*Amp*f_coeffTEx_z[pos,k]*css_z
            #ETEy_z[i,j]+=H0*Amp*f_coeffTEy_z[pos,k]*scs_z

            #TM modes x
            #ETMy_x[i,j]+=E0*Amp*f_coeffTMy_x[pos,k]*css_x
            #ETMz_x[i,j]+=E0*Amp*f_coeffTMz_x[pos,k]*scs_x
            ETMx_x[i,j]+=E0*Amp*f_coeffTMx_x[pos,k]*ssc_x

            #TE modes x
            #ETEy_x[i,j]+=H0*Amp*f_coeffTEy_x[pos,k]*css_x
            #ETEz_x[i,j]+=H0*Amp*f_coeffTEz_x[pos,k]*scs_x

            ##TM modes y
            #ETMz_y[i,j]+=E0*Amp*f_coeffTMz_y[pos,k]*css_y
            ETMx_y[i,j]+=E0*Amp*f_coeffTMx_y[pos,k]*scs_y
            #ETMy_y[i,j]+=E0*Amp*f_coeffTMy_y[pos,k]*ssc_y

            ##TE modes y
            #ETEz_y[i,j]+=H0*Amp*f_coeffTEz_y[pos,k]*css_y
            ETEx_y[i,j]+=H0*Amp*f_coeffTEx_y[pos,k]*scs_y
          end

          #if M[k]>0 && N[k]>0 && P[k]==0
            #ETMz_z[i,j] += E0*Amp*f_coeffTMz_z[pos,k]*ssc_z #TM mode z
            #ETEz_y[i,j]+=H0*Amp*f_coeffTEz_y[pos,k]*css_y #TE mode y
            #ETEz_x[i,j]+=H0*Amp*f_coeffTEz_x[pos,k]*scs_x #TE mode x
          #end

          if M[k]==0 && N[k]>0 && P[k]>0
            ETEx_z[i,j] += H0*Amp*f_coeffTEx_z[pos,k]*css_z #TE mode z
            ETMx_x[i,j]+=E0*Amp*f_coeffTMx_x[pos,k]*ssc_x #TM mode x
            ETEx_y[i,j]+=H0*Amp*f_coeffTEx_y[pos,k]*scs_y #TE mode y
          end

          #if M[k]>0 && N[k]==0 && P[k]>0
            #ETEy_z[i,j]+=H0*Amp*f_coeffTEy_z[pos,k]*scs_z #TE mode z
            #ETMy_y[i,j]+=E0*Amp*f_coeffTMy_y[pos,k]*ssc_y #TM mode y
            #ETEy_x[i,j]+=H0*Amp*f_coeffTEy_x[pos,k]*css_x #TE mode x
          #end
        end
      end
    end
    Ex[u,pos,:]=ETMx_x+ETMx_y+ETMx_z+ETEx_y+ETEx_z
    #Ey[u,pos,:]=ETMy_x+ETMy_y+ETMy_z+ETEy_x+ETEy_z
    #Ez[u,pos,:]=ETMz_x+ETMz_y+ETMz_z+ETEz_y+ETEz_x
  end
  M=[]
  N=[]
  P=[]
end
u=1
cmax=maximum(abs(Ex[u,:,:,:]))

for o=1:n_steps
  figure(num=2,figsize=(8,5))
  ftext=round(f[u]/1e6*100)/100
  xtext1=round(x_b[o]*10)/10
  ytext1=round(y_b[o]*10)/10
  ztext1=round(z_b[o]*10)/10
  #title("Emitter at \$($xtext1, $ytext1, $ztext1)\$",fontsize=17)
  contourf(z,y,squeeze(squeeze(abs(Ex[u,o,:,:]),1),1),10,cmap="RdBu_r")
  plot(z_b[o],y_b[o],marker="o",markerfacecolor="w", markersize=10)
  clim(0,cmax)
  if o==1
    colorbar(format="%.2f",orientation="horizontal")
  end
  xlabel("\$ z\$/m",fontsize=10)
  ylabel("\$ y\$/m",fontsize=10)
  axis("image")
  grid()
  savefig("Exk_$o.pdf",bbox_inches="tight")
  clf()
end



