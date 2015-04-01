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
const f0=40e6
const f1=1000e6

#rf=.01
#nf=int(1/log(1+rf)*log(f1/f0))
#f=zeros(nf)
#for i=1:nf
#  f[i]=f0*(1+rf)^i
#end

const df=1e6
const nf=1#int((f1-f0)/df)+1
const f=511e6#linspace(f0,f1,nf)#array([20e6,40e6,45e6,56e6,65e6])#linspace(20e6,.3e9,281)#linspace(20e6,3e8,1121)
const V=a*b*d
const Q4=16*pi^2*V*f.^3/(1*c^3)
const Q1=3*V./(2*1*(2./sqrt(2*pi*f.*mur*mu0*1*1e5)))/200
const Q=1./(1./Q1+1./Q4)

const x=1
dy=0.02
dz=0.02
ny=int(b/dy)+1
nz=int(d/dz)+1
const y=linspace(0,b,ny)
const z=linspace(0,d,nz)


n_steps=360
phi=linspace(2*pi/n_steps,2*pi,n_steps)
n_sources = 1

x_b=zeros(n_sources,n_steps)
y_b=zeros(n_sources,n_steps)
z_b=zeros(n_sources,n_steps)

for i=1:n_sources
  phi0=rand()*2pi
  r0=1#rand()
  x_b[i,:]=ones(n_steps)*1#rand()*(a-1)+0.5
  y_b[i,:]=1.5+r0*sin(phi+phi0)
  z_b[i,:]=1.5+r0*cos(phi+phi0)
end


E0=1
H0=1/(c*mu0*mur)

Nmodesx=zeros(nf)
Nx=zeros(nf)



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
  println("$(f[u]), $nm resonant freq.")
  f_coeffTMx_x=zeros(n_steps,nm)
  f_coeffTMx_y=zeros(n_steps,nm)
  f_coeffTMx_z=zeros(n_steps,nm)
  f_coeffTEx_y=zeros(Complex128,n_steps,nm)
  f_coeffTEx_z=zeros(Complex128,n_steps,nm)

  f_coeffTMy_x=zeros(n_steps,nm)
  f_coeffTMy_y=zeros(n_steps,nm)
  f_coeffTMy_z=zeros(n_steps,nm)
  f_coeffTEy_x=zeros(Complex128,n_steps,nm)
  f_coeffTEy_z=zeros(Complex128,n_steps,nm)

  f_coeffTMz_x=zeros(n_steps,nm)
  f_coeffTMz_y=zeros(n_steps,nm)
  f_coeffTMz_z=zeros(n_steps,nm)
  f_coeffTEz_y=zeros(Complex128,n_steps,nm)
  f_coeffTEz_x=zeros(Complex128,n_steps,nm)
  for pos=1:n_steps
    Ex=zeros(Complex128,ny,nz)
    Ey=zeros(Complex128,ny,nz)
    Ez=zeros(Complex128,ny,nz)
    ETMx_x=zeros(ny,nz)
    ETMx_y=zeros(ny,nz)
    ETMx_z=zeros(ny,nz)
    ETEx_y=zeros(Complex128,ny,nz)
    ETEx_z=zeros(Complex128,ny,nz)

    ETMy_x=zeros(ny,nz)
    ETMy_y=zeros(ny,nz)
    ETMy_z=zeros(ny,nz)
    ETEy_x=zeros(Complex128,ny,nz)
    ETEy_z=zeros(Complex128,ny,nz)

    ETMz_x=zeros(ny,nz)
    ETMz_y=zeros(ny,nz)
    ETMz_z=zeros(ny,nz)
    ETEz_y=zeros(Complex128,ny,nz)
    ETEz_x=zeros(Complex128,ny,nz)

    xe=x_b[:,pos]
    ye=y_b[:,pos]
    ze=z_b[:,pos]
    #println(" Step $pos/$n_steps")
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
        f_coeffTMx_z[pos,k] = -E0*Amp*kx*kz/(k2mnp-kz^2)*css_ze
        f_coeffTMy_z[pos,k] = E0*Amp*ky*kz/(k2mnp-kz^2)*scs_ze
        f_coeffTMz_z[pos,k] = E0*Amp*ssc_ze
        #TE modes z
        f_coeffTEx_z[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze
        f_coeffTEy_z[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze

        #TM modes x
        f_coeffTMy_x[pos,k] = -E0*Amp*ky*kx/(k2mnp-kx^2)*css_xe
        f_coeffTMz_x[pos,k] = E0*Amp*kz*kx/(k2mnp-kx^2)*scs_xe
        f_coeffTMx_x[pos,k] = E0*Amp*ssc_xe
        #TE modes x
        f_coeffTEy_x[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe
        f_coeffTEz_x[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe

        ##TM modes y
        f_coeffTMz_y[pos,k] += -H0*Amp*E0*Amp*kz*ky/(k2mnp-ky^2)*css_ye
        f_coeffTMx_y[pos,k] = H0*Amp*E0*Amp*kx*ky/(k2mnp-ky^2)*scs_ye
        f_coeffTMy_y[pos,k] = H0*Amp*E0*Amp*ssc_ye
        ##TE modes y
        f_coeffTEz_y[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye
        f_coeffTEx_y[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye
      end


      if M[k]>0 && N[k]>0 && P[k]==0
        f_coeffTMz_z[pos,k] = E0*Amp*ssc_ze #TM mode z
        f_coeffTEz_y[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye #TE mode y
        f_coeffTEz_x[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe #TE mode x
      end

      if M[k]==0 && N[k]>0 && P[k]>0
        f_coeffTEx_z[pos,k] += -H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze #TE mode z
        f_coeffTMx_x[pos,k] += E0*Amp*ssc_xe #TM mode x
        f_coeffTEx_y[pos,k] += H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye #TE mode y
      end

      if M[k]>0 && N[k]==0 && P[k]>0
        f_coeffTEy_z[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze #TE mode z
        f_coeffTMy_y[pos,k] = E0*Amp*ssc_ye #TM mode y
        f_coeffTEy_x[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe #TE mode x
      end

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
            ETMx_z[i,j] += f_coeffTMx_z[pos,k]*css_z
            ETMy_z[i,j] += f_coeffTMy_z[pos,k]*scs_z
            ETMz_z[i,j] += f_coeffTMz_z[pos,k]*ssc_z

            #TE modes z
            ETEx_z[i,j] += f_coeffTEx_z[pos,k]*css_z
            ETEy_z[i,j] += f_coeffTEy_z[pos,k]*scs_z

            #TM modes x
            ETMy_x[i,j] += f_coeffTMy_x[pos,k]*css_x
            ETMz_x[i,j] += f_coeffTMz_x[pos,k]*scs_x
            ETMx_x[i,j] += f_coeffTMx_x[pos,k]*ssc_x

            #TE modes x
            ETEy_x[i,j] += f_coeffTEy_x[pos,k]*css_x
            ETEz_x[i,j] += f_coeffTEz_x[pos,k]*scs_x

            #TM modes y
            ETMz_y[i,j] += f_coeffTMz_y[pos,k]*css_y
            ETMx_y[i,j] += f_coeffTMx_y[pos,k]*scs_y
            ETMy_y[i,j] += f_coeffTMy_y[pos,k]*ssc_y

            #TE modes y
            ETEz_y[i,j] += f_coeffTEz_y[pos,k]*css_y
            ETEx_y[i,j] += f_coeffTEx_y[pos,k]*scs_y
          end

          if M[k]>0 && N[k]>0 && P[k]==0
            ETMz_z[i,j] += f_coeffTMz_z[pos,k]*ssc_z #TM mode z
            ETEz_y[i,j] += f_coeffTEz_y[pos,k]*css_y #TE mode y
            ETEz_x[i,j] += f_coeffTEz_x[pos,k]*scs_x #TE mode x
          end

          if M[k]==0 && N[k]>0 && P[k]>0
            ETEx_z[i,j] += f_coeffTEx_z[pos,k]*css_z #TE mode z
            ETMx_x[i,j] += f_coeffTMx_x[pos,k]*ssc_x #TM mode x
            ETEx_y[i,j] += f_coeffTEx_y[pos,k]*scs_y #TE mode y
          end

          if M[k]>0 && N[k]==0 && P[k]>0
            ETEy_z[i,j] += f_coeffTEy_z[pos,k]*scs_z #TE mode z
            ETMy_y[i,j] += f_coeffTMy_y[pos,k]*ssc_y #TM mode y
            ETEy_x[i,j] += f_coeffTEy_x[pos,k]*css_x #TE mode x
          end
        end
      end
    end
    Ex=ETMx_x+ETMx_y+ETMx_z+ETEx_y+ETEx_z
    Ey=ETMy_x+ETMy_y+ETMy_z+ETEy_x+ETEy_z
    Ez=ETMz_x+ETMz_y+ETMz_z+ETEz_y+ETEz_x

    figure(num=2,figsize=(21.12,11.52))
    ftext=round(f[u]/1e6*100)/100
    xtext1=round(x_b[pos]*10)/10
    ytext1=round(y_b[pos]*10)/10
    ztext1=round(z_b[pos]*10)/10
    subplot(131)
    title("\$ \\vert E_x\\vert\$",fontsize=17)
    pcolor(y,z,abs(Ex'),cmap="RdBu_r")
    plot(y_b[pos],z_b[pos],marker="o",markerfacecolor="w", markersize=10)
    clim(0,0.3)
    #colorbar(format="%.2f")
    xlabel("\$ y\$/m",fontsize=15)
    ylabel("\$ z\$/m",fontsize=15)
    axis("image")
    grid()

    subplot(132)
    title("\$ \\vert E_y\\vert\$",fontsize=17)
    pcolor(y,z,abs(Ey'),cmap="RdBu_r")
    plot(y_b[pos],z_b[pos],marker="o",markerfacecolor="w", markersize=10)
    clim(0,0.3)
    #colorbar(format="%.2f")
    xlabel("\$ y\$/m",fontsize=15)
    #ylabel("\$ z\$/m",fontsize=15)
    axis("image")
    grid()

    subplot(133)
    title("\$ \\vert E_z\\vert\$",fontsize=17)
    pcolor(y,z,abs(Ez'),cmap="RdBu_r")
    plot(y_b[pos],z_b[pos],marker="o",markerfacecolor="w", markersize=10)
    clim(0,0.3)
    colorbar(format="%.2f")
    xlabel("\$ y\$/m",fontsize=15)
    #ylabel("\$ z\$/m",fontsize=15)
    axis("image")
    grid()
    savefig("Exyz_$pos.png")
    clf()

end

end
