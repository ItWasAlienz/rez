/**
	@file
	rez~ - a looper thang

	@ingroup	MSP
*/

#include "ext.h"
#include "z_dsp.h"
#include "math.h"
#include "ext_buffer.h"
#include "ext_atomic.h"
#include "ext_obex.h"
#include "ext_obex_util.h"

#ifdef MAC_VERSION
#define FRCNLN inline __attribute__((always_inline))
#endif
#ifdef WIN_VERSION
#define FRCNLN __forceinline
#endif

typedef struct _rez {
	t_pxobject obj;
	t_buffer_ref *buf;
	t_symbol *bufname;
    t_symbol *mode;
    t_ptr_int chan;
    t_ptr_int bchan;
    t_ptr_int bframes;
    t_ptr_int pstart;
    t_ptr_int pend;
    t_ptr_int rstart;
    t_ptr_int rend;
    t_ptr_int nstart;
    t_ptr_int nend;
    t_ptr_int nrstart;
    t_ptr_int nrend;
    t_ptr_int fade;
    t_ptr_int pfad;
    t_ptr_int rfad;
    t_ptr_int rxfad;
    t_ptr_int xfad;
    t_ptr_int rspprev;
    t_double rphprev;
    t_double dpos;
    t_double xpos;
    t_double rpos;
    t_double rxpos;
    t_double phprev;
    t_double spprev;
    t_double sr;
    t_double srscale;
    t_bool wr;
    t_bool rwr;
    t_bool pparamchange;
    t_bool rparamchange;
    t_bool inOderAus;
    t_bool rinOderAus;
	t_bool buf_modd;
    t_bool wram;
    t_bool oneshot;
    t_bool nshot;
    t_bool rnshot;
} t_rez;


void *rez_new(t_symbol *s,  long argc, t_atom *argv);
void rez_free(t_rez *x);
t_max_err rez_notify(t_rez *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
t_max_err rez_attrfad_get(t_rez *x, t_object *attr, long *argc, t_atom **argv);
t_max_err rez_attrfad_set(t_rez *x, t_object *attr, long argc, t_atom *argv);
t_max_err rez_attr1shot_get(t_rez *x, t_object *attr, long *argc, t_atom **argv);
t_max_err rez_attr1shot_set(t_rez *x, t_object *attr, long argc, t_atom *argv);
void rez_assist(t_rez *x, void *b, long m, long a, char *s);
void rez_fade(t_rez *x, t_ptr_int fd);
void rez_limits(t_rez *x);
void rez_set(t_rez *x, t_symbol *s, long ac, t_atom *av);
void rez_doset(t_rez *x, t_symbol *s, long ac, t_atom *av);
void rez_pstart(t_rez *x, t_double ps);
void rez_pend(t_rez *x, t_double pe);
void rez_rstart(t_rez *x, t_double rs);
void rez_rend(t_rez *x, t_double re);
void rez_dblclick(t_rez *x);
void rez_sperform64(t_rez *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void rez_psperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp);
void rez_rsperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp);
void rez_mperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp);
void rez_pmperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp);
void rez_rmperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp);
void rez_dsp64(t_rez *x, t_object *d64, short *count, double sr, long maxvecsize, long flgs);

static t_class *s_rez_class;
static t_symbol *ps_buffer_modified, *ps_mode_a, *ps_mode_b, *ps_nada;

//easing function(fade/crossfade)
static FRCNLN double eas_func_up(double y1, double ramp, long fad)
{ return y1*(0.5*(1.0-cos((1.0-(((double)fad)/ramp))*PI))); }

static FRCNLN double eas_func_dwn(double y1, double ramp, long fad)
{ return y1*(0.5*(1.0-cos((((double)fad)/ramp)*PI))); }

static FRCNLN double eas_rec_up(double y1, double rmp, double odb, long fad)
{ return y1*(1.0-((odb*0.5)*(1.0-cos((1.0-((double)fad/rmp))*PI)))); }

static FRCNLN double eas_rec_dwn(double y1, double rmp, double odb, long fad)
{ return y1*(1.0-((odb*0.5)*(1.0-cos(((double)fad/rmp)*PI)))); }

static FRCNLN void interp_index(t_ptr_int p, t_ptr_int *ndxz,t_ptr_int *ndxb,t_ptr_int *ndxc, t_double spd,t_ptr_int frms)
{
    t_ptr_int dr = (spd>=0)?1:-1;
    *ndxz = p - dr; frms -= 1;                              //index-calc for cubic interpolation
    if(*ndxz < 0) *ndxz = frms + *ndxz; else if(*ndxz > frms) *ndxz = *ndxz - frms;

    *ndxb = p + dr;
    if(*ndxb < 0) *ndxb = frms + *ndxb; else if(*ndxb > frms) *ndxb = *ndxb - frms;

    *ndxc = *ndxb + dr;
    if(*ndxc < 0) *ndxc = frms + *ndxc; else if(*ndxc > frms) *ndxc = *ndxc - frms;
    return;
}

static FRCNLN double interp(double f, double z, double a, double b, double c)
{ return (!f) ? a : ((((0.5*(c-z) + 1.5*(a-b))*f + (z - 2.5*a + b + b - 0.5*c))*f + (0.5*(b-z)))*f + a); }

void ext_main(void *r)
{
	t_class *c = class_new("rez~", (method)rez_new, (method)rez_free, sizeof(t_rez), NULL, A_GIMME, 0);

	class_addmethod(c, (method)rez_dsp64,	"dsp64",	A_CANT, 0);
    class_addmethod(c, (method)rez_assist,  "assist",   A_CANT, 0);
    class_addmethod(c, (method)rez_dblclick,"dblclick", A_CANT, 0);
    class_addmethod(c, (method)rez_notify,  "notify",   A_CANT, 0);
	class_addmethod(c, (method)rez_pstart,	"in",	A_FLOAT, 0);
	class_addmethod(c, (method)rez_pend,	"aus",		A_FLOAT, 0);
    class_addmethod(c, (method)rez_rstart,  "rin",   A_FLOAT, 0);
    class_addmethod(c, (method)rez_rend,    "raus",     A_FLOAT, 0);
	class_addmethod(c, (method)rez_set,		"set",		A_GIMME, 0);

    CLASS_ATTR_LONG(c, "fade", 0, t_rez, fade);
    CLASS_ATTR_SAVE(c, "fade", 0);
    CLASS_ATTR_ACCESSORS(c, "fade", (method)rez_attrfad_get, (method)rez_attrfad_set);
    
    CLASS_ATTR_LONG(c, "oneshot", 0, t_rez, oneshot);
    CLASS_ATTR_SAVE(c, "oneshot", 0);
    CLASS_ATTR_ACCESSORS(c, "oneshot", (method)rez_attr1shot_get, (method)rez_attr1shot_set);
    
	class_dspinit(c); class_register(CLASS_BOX, c); s_rez_class = c;

	ps_buffer_modified = gensym("buffer_modified"); ps_nada = gensym("");
    ps_mode_a = gensym("play"); ps_mode_b = gensym("rec");
}

void *rez_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_rez *x = (t_rez *)object_alloc(s_rez_class);
    if(x)
    {
        t_symbol *bufname=0; t_symbol *mode=0; t_ptr_int chan=2;
        bufname = atom_getsymarg(0,argc,argv);
        chan = atom_getintarg(1,argc,argv); if(!chan)chan=2;
        mode = atom_getsymarg(2,argc,argv);
    
        if(mode == ps_mode_a){ dsp_setup((t_pxobject *)x,2); }else if(mode == ps_mode_b)
        { if(chan>1) dsp_setup((t_pxobject *)x,5); else dsp_setup((t_pxobject *)x,4); }
        else{ if(chan>1) dsp_setup((t_pxobject *)x,7); else dsp_setup((t_pxobject *)x,6); }
    
        x->mode=mode; x->bufname=bufname; x->chan=chan; x->fade=128;
        x->rparamchange=x->pparamchange=x->rinOderAus=x->inOderAus=1;
        x->rxfad=x->xfad=x->pfad=x->rfad=-1; x->wram=x->rspprev=x->buf_modd=x->rnshot=0;
        x->nshot=x->nstart=x->nrstart=x->rstart=x->pstart=x->rwr=x->wr=x->oneshot=0;
        x->rphprev=x->rpos=x->rxpos=x->xpos=x->dpos=x->spprev=x->phprev=0.;
    
        if(mode!=ps_mode_b)
        {
            if(chan>1)outlet_new((t_object *)x, "signal");     // audio outlet
            outlet_new((t_object *)x, "signal");        // audio outlet
        }
        if(mode==ps_nada){ outlet_new((t_object *)x, "signal"); }
        outlet_new((t_object *)x, "signal");       // phase outlet

        // create new buffer reference, initially a buffer with provided-name
        x->buf = buffer_ref_new((t_object *)x, bufname); x->obj.z_misc |= Z_NO_INPLACE;
        attr_args_process(x, argc, argv);
    }
    return (x);
}
// add object to dsp-chain(and call appropriate perform method; setup according to buffer~)
void rez_dsp64(t_rez *x, t_object *d64, short *count, double sr, long maxvecsize, long flgs)
{
    x->sr = sr;
    t_atom iav[2];
    atom_setsym(iav, x->bufname);
    atom_setlong(iav+1, x->chan);
    rez_set(x, x->bufname, 2, iav);
    if(x->chan>1)
    {
        if(x->mode==ps_mode_a)
            object_method(d64, gensym("dsp_add64"), x, rez_psperform64, 0, NULL);
        else if(x->mode==ps_mode_b)
            object_method(d64, gensym("dsp_add64"), x, rez_rsperform64, 0, NULL);
        else
            object_method(d64, gensym("dsp_add64"), x, rez_sperform64, 0, NULL);
    }
    else
    {
        if(x->mode==ps_mode_a)
            object_method(d64, gensym("dsp_add64"), x, rez_pmperform64, 0, NULL);
        else if(x->mode==ps_mode_b)
            object_method(d64, gensym("dsp_add64"), x, rez_rmperform64, 0, NULL);
        else
            object_method(d64, gensym("dsp_add64"), x, rez_mperform64, 0, NULL);
    }
}
// free buffer~ memory, and object from dsp-chain
void rez_free(t_rez *x){ dsp_free((t_pxobject *)x);  object_free(x->buf); }
// handles necessary notifications for when buffer appears/disappears/is-modded.
t_max_err rez_notify(t_rez *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == ps_buffer_modified) x->buf_modd = true;
	return buffer_ref_notify(x->buf, s, msg, sender, data);
}

void rez_assist(t_rez *x, void *b, long m, long a, char *s)
{
    t_symbol *mode = x->mode; t_ptr_int chan= x->chan;
	if (m == ASSIST_INLET) // inlets
    {
        if(mode == ps_mode_a)
        {
            switch (a)
            {
                case 0: snprintf_zero(s,256,"(signal)starting-phase(from 0. to 1.)"); break;
                case 1: snprintf_zero(s, 256,"(signal)playback speed(from -∞. to +∞.)"); break;
            }
        }
        else if(mode == ps_mode_b)
        {
            if(chan>1)
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)starting-rec-phase(from 0. to 1.)"); break;
                    case 1: snprintf_zero(s, 256,"(signal)recording-speed(-1.,0.,or +1. only)"); break;
                    case 2: snprintf_zero(s, 256, "(signal)overdub-amount(0. to 1.)"); break;
                    case 3: snprintf_zero(s, 256, "(signal)left-recording-input"); break;
                    case 4: snprintf_zero(s, 256, "(signal)right-recording-input"); break;
                }
            }
            else
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)starting-rec-phase(from 0. to 1.)"); break;
                    case 1: snprintf_zero(s, 256,"(signal)recording-speed(-1.,0.,or +1. only)"); break;
                    case 2: snprintf_zero(s, 256, "(signal)overdub-amount(0. to 1.)"); break;
                    case 3: snprintf_zero(s, 256, "(signal)left-recording-input"); break;
                }
            }
        }
        else
        {
            if(chan>1)
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)starting-phase(from 0. to 1.)"); break;
                    case 1: snprintf_zero(s, 256,"(signal)playback speed(from -∞. to +∞.)"); break;
                    case 2: snprintf_zero(s,256,"(signal)starting-rec-phase(from 0. to 1.)"); break;
                    case 3: snprintf_zero(s, 256,"(signal)recording-speed(-1.,0.,or +1. only)"); break;
                    case 4: snprintf_zero(s, 256,"(signal)overdub-amount(0. to 1.)"); break;
                    case 5: snprintf_zero(s, 256,"(signal)left-recording-input"); break;
                    case 6: snprintf_zero(s, 256,"(signal)right-recording-input"); break;
                }
            }
            else
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)starting-phase(from 0. to 1.)"); break;
                    case 1: snprintf_zero(s, 256,"(signal)playback speed(from -∞. to +∞.)"); break;
                    case 2: snprintf_zero(s,256,"(signal)starting-rec-phase(from 0. to 1.)"); break;
                    case 3: snprintf_zero(s, 256,"(signal)recording-speed(-1.,0.,or +1. only)"); break;
                    case 4: snprintf_zero(s, 256,"(signal)overdub-amount(0. to 1.)"); break;
                    case 5: snprintf_zero(s, 256,"(signal)left-recording-input"); break;
                }
            }
        }
	}
	else
    {
        if(mode == ps_mode_a)
        {
            if(chan>1)
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)left-output"); break;
                    case 1: snprintf_zero(s, 256,"(signal)right-output"); break;
                    case 2: snprintf_zero(s, 256,"(signal)playback-phase(0. to 1.)"); break;
                }
            }
            else
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)left-output"); break;
                    case 1: snprintf_zero(s, 256,"(signal)playback-phase(0. to 1.)"); break;
                }
            }
        }
        else if(mode == ps_nada)
        {
            if(chan>1)
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)left-output"); break;
                    case 1: snprintf_zero(s, 256,"(signal)right-output"); break;
                    case 2: snprintf_zero(s, 256,"(signal)playback-phase(0. to 1.)"); break;
                    case 3: snprintf_zero(s, 256,"(signal)recording-phase(0. to 1.)"); break;
                }
            }
            else
            {
                switch (a)
                {
                    case 0: snprintf_zero(s,256,"(signal)left-output"); break;
                    case 1: snprintf_zero(s, 256,"(signal)playback-phase(0. to 1.)"); break;
                    case 2: snprintf_zero(s, 256,"(signal)recording-phase(0. to 1.)"); break;
                }
            }
        }else{ switch(a){ case 0: snprintf_zero(s, 256,"(signal)recording-phase(0. to 1.)"); } }
    }
}

t_max_err rez_attrfad_set(t_rez *x, t_object *attr, long argc, t_atom *argv)
{ t_ptr_int fade = atom_getlong(argv); x->fade = (fade<0)?0:((fade>65536)?65536:fade); return 0; }

t_max_err rez_attrfad_get(t_rez *x, t_object *attr, long *argc, t_atom **argv)
{
    char alloc; t_ptr_int fade = x->fade;
    atom_alloc(argc, argv, &alloc); atom_setlong(*argv, fade); return 0;
}

t_max_err rez_attr1shot_set(t_rez *x, t_object *attr, long argc, t_atom *argv)
{ t_bool oneshot = atom_getlong(argv); x->oneshot = oneshot; return 0; }

t_max_err rez_attr1shot_get(t_rez *x, t_object *attr, long *argc, t_atom **argv)
{
    char alloc; t_bool oneshot = x->oneshot;
    atom_alloc(argc, argv, &alloc); atom_setlong(*argv, oneshot); return 0;
}

void rez_limits(t_rez *x)
{
    t_double sr = x->sr;
	t_buffer_obj *b = buffer_ref_getobject(x->buf); // get the actual buffer object from our ref
	if (b)
    {
		x->bchan=buffer_getchannelcount(b); // buffer~'s number of chans
        x->srscale=buffer_getsamplerate(b)/sr;
        x->bframes=buffer_getframecount(b); // buffer~-length of 1 channel
    }
}

void rez_doset(t_rez *x, t_symbol *s, long ac, t_atom *av)
{
	t_symbol *name; name = (ac) ? atom_getsym(av) : gensym("");
	buffer_ref_set(x->buf, name); rez_limits(x);
}

// calls setting buffer ref be on main thread only
void rez_set(t_rez *x,t_symbol *s,long ac,t_atom *av){ defer(x,(method)rez_doset,s,ac,av); }

//___________LOOP_BOUNDARY_SETUP FOR PLAYBACK AND RECORDING___________
//.."wr"/"rwr" = wraparound flag is on if "in"(pstart)/"rin"(rstart) > "aus"(pend)/"raus"(rend)...
//..With wraparound, playback plays through end-to-start of buffer~ to reach "aus"(pend)/"raus"(rend)
//_____"pparamchange"/"rparamchange" are flags to tell the perform method to register changes...
//.....(once safely within distance for proper observance of loop-boundaries/crossfades/fades,..
//..perform-meth writes to "pparamchange"/"rparamchange" tracking internal registration of changes)
void rez_pstart(t_rez *x, t_double ps)
{
    ps = (ps>1.0)?1.0:((ps<0.0)?0.0:ps); x->pstart = ps*(x->bframes-1);
    if(x->pstart > x->pend) x->wr=1; else x->wr=0; x->pparamchange = 1;
}

void rez_pend(t_rez *x, t_double pe)
{
    pe = (pe>1.0)?1.0:((pe<0.0)?0.0:pe); x->pend = pe*(x->bframes-1);
    if(x->pstart > x->pend) x->wr=1; else x->wr=0; x->pparamchange = 1;
}

void rez_rstart(t_rez *x, t_double rs)
{
    rs = (rs>1.0)?1.0:((rs<0.0)?0.0:rs); x->rstart = rs*(x->bframes-1);
    if(x->rstart > x->rend) x->rwr=1; else x->rwr=0; x->rparamchange = 1;
}

void rez_rend(t_rez *x, t_double re)
{
    re = (re>1.0)?1.0:((re<0.0)?0.0:re); x->rend = re*(x->bframes-1);
    if(x->rstart > x->rend) x->rwr=1; else x->rwr=0; x->rparamchange = 1;
}

void rez_dblclick(t_rez *x) { buffer_view(buffer_ref_getobject(x->buf)); }

void rez_sperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
	t_double		*ph = is[0];
    t_double        *spd = is[1];
    t_double        *rph = is[2];
    t_double        *rspd = is[3];
    t_double        *ovdb = is[4];
    t_double        *rinL = is[5];
    t_double        *rinR = is[6];
	t_double		*outL = os[0];
    t_double        *outR = os[1];
    t_double        *outPh = os[2];
    t_double        *outRPh = os[3];

	long		 n=smps;                                           float *b;
    t_double	 phprev,inpos,dpos,sp,rcL,rcR,frac,dfrnc,rpos,rxpos,rsp;
    t_double     spprev,oL,oR,oPh,oRPh,xL,xR,xpos,scldsr,odb,rphprev,rinpos;
    t_ptr_int    rndx,rxdx,bnc,bframes,pstart,pend,rstart,rend,nstart,nend,nrstart,nrend;
    t_ptr_int    indx,indxz,indxb,indxc,xfad,pfad,rfad,rxfad,fade,rs,rspprev;
    t_bool       wr,rwr,pparamch,rparamch,inOderAus,rinOderAus,wram,oneshot,nshot,rnshot;
	t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

	b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot; nshot=x->nshot;
    bnc=x->bchan; bframes=x->bframes-1; pstart=x->pstart; pend=x->pend; inOderAus=x->inOderAus;
    pparamch=x->pparamchange; pfad=x->pfad; xfad=x->xfad; fade=x->fade; xpos=x->xpos;
    rpos=x->rpos; wr=x->wr;  rparamch=x->rparamchange; rfad=x->rfad;
    rstart=x->rstart; rend=x->rend; rnshot=x->rnshot; phprev=x->phprev; rphprev=x->rphprev;
    spprev=x->spprev;  dpos=x->dpos; wram=x->wram; scldsr=x->srscale;
    nend=x->nend; nstart=x->nstart; nrend=x->nrend; nrstart=x->nrstart;
    rxfad=x->rxfad; rxpos=x->rxpos; rspprev=x->rspprev; rwr=x->rwr; rinOderAus=x->rinOderAus;
    
    while(n--)
    {
        inpos=*ph++; sp=*spd++; rinpos=*rph++; rsp=*rspd++; odb=*ovdb++; rcL=*rinL++; rcR=*rinR++;
        inpos=(inpos<0.0)?0.0:((inpos>1.0)?1.0:inpos); inpos*=bframes; sp*=scldsr;
        rinpos=(rinpos<0.0)?0.0:((rinpos>1.0)?1.0:rinpos); rinpos*=bframes;
        if(rsp>0.0) rs=1; else if(rsp<0.0) rs=-1; else rs=0;
//inOderAus(from SachaBaronCohen's "Bruno" skit)-flag indicates if position is 'in'/'out' of a loop..
//..if 'in' loop-boundary,save messages(pstart/pend/rstart/rend) internally(nstart/nend/nrstart/nrend)
        if(oneshot)//ONESHOT_SECTION
        {               //..........................ONESHOT-RECORDING
            if((rspprev==0)&&(rs!=0))//on speed-change for 'on',check position from loop bounds..
            {       //..setup fade-in/fade-out,internal one-shot flag(nshot),& pos/speed changes
                rpos=rinpos; rfad=fade-1; wram=1; rnshot=1; rspprev=rs;
                if(rwr)
                {
                    if((rpos>(rend-fade))&&(rpos<(rstart+fade))){ rinOderAus=0; }
                    else { rinOderAus=1; nrstart=rstart; nrend=rend; }
                }
                else
                {
                    if((rpos>(rend-fade))||(rpos<(rstart+fade))){ rinOderAus=0; }
                    else { rinOderAus=1; nrstart=rstart; nrend=rend; }
                }
            }
            else if(((rspprev!=0)&&(rs==0)) && rnshot){ rfad=fade-1; rnshot=0; }
            
            if(rnshot)//rnshot flag determines if the state is fading-in-recording/recording-regular..
            {
                if(rwr)
                {
                    if((rpos<(rend-fade))||(rpos>(rstart+fade)))
                    { rinOderAus=1; nrstart=rstart; nrend=rend; }
                }
                else
                {
                    if((rpos<(rend-fade))&&(rpos>(rstart+fade)))
                    { rinOderAus=1; nrstart=rstart; nrend=rend; }
                }
                if(rinOderAus)
                {
                    if(rwr)
                    { if((rpos>=(nrend-fade))&&(rpos<=(nrstart+fade))){ rfad=fade-1; rnshot=0; } }
                    else
                    { if((rpos>=(nrend-fade))||(rpos<=(nrstart+fade))){ rfad=fade-1; rnshot=0; } }
                }
                        //buffer~-boundary constraint, and positional-increment
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=rpos; rpos+=rs;
                if(rfad>=0)//Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--;
                }else{ b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); b[rndx*bnc+1]=rcR+(b[rndx*bnc+1]*odb); };
            }
            else    //..or(when rnshot is 0, the state is) fading-out-recording/stopped-recording
            {
                if(rfad>=0)
                {                 //....buffer~-boundary constraint, and positional-increment
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=rpos; rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_dwn(rcR,fade,rfad)+eas_rec_dwn(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; rspprev=rs; }
                } else { rspprev=rs; }
            }
            
            //......................................ONESHOT-PLAYBACK
            if((spprev==0.0)&&(sp!=0.0))//on speed-change for 'on',check position from loop bounds..
            {     //..setup fade-in/fade-out,internal one-shot flag(nshot),& position/speed changes
                dpos=inpos; pfad=fade-1; nshot=1; spprev=sp;
                if(wr)
                {
                    if((dpos>(pend-fade))&&(dpos<(pstart+fade))){ inOderAus=0; }
                    else { inOderAus=1; nstart=pstart; nend=pend; }
                }
                else{ if((dpos>(pend-fade))||(dpos<(pstart+fade)))
                { inOderAus=0; }else{ inOderAus=1; nstart=pstart; nend=pend; } }
            }
            else if(((spprev!=0.0)&&(sp==0.0)) && nshot){ pfad=fade-1; nshot=0; }
            
            if(nshot)    //nshot flag determines if the state is fading-in/playing-back-regular..
            {
                if(wr)
                {
                    if((dpos<(pend-fade))||(dpos>(pstart+fade)))
                    { inOderAus=1; nstart=pstart; nend=pend; }
                }
                else
                {
                    if((dpos<(pend-fade))&&(dpos>(pstart+fade)))
                    { inOderAus=1; nstart=pstart; nend=pend; }
                }
                if(inOderAus)
                {
                    if(wr){ if((dpos>=(nend-fade))&&(dpos<=(nstart+fade))){ pfad=fade-1; nshot=0; } }
                    else{ if((dpos>=(nend-fade))||(dpos<=(nstart+fade))){ pfad=fade-1; nshot=0; } }
                }
                //buffer~-boundary constraint, positional-increment, and bicubic interpolation
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                if(sp>0){frac=dpos-indx;}else if(sp<0){frac=1.0-(dpos-indx);}else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                //fade-in
                if(pfad>=0){ oL=eas_func_up(xL,fade,pfad); oR=eas_func_up(xR,fade,pfad); pfad--; }
                else{ oL=xL; oR=xR; }
            }
            else        //..or(when nshot is 0, the state is) fading-out/silenced
            {//fade-Out
                if(pfad>=0)
                { //buffer~-boundary constraint, positional-increment, and bicubic interpolation
                    if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                    if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                    dpos+=spprev;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_dwn(xL, fade, pfad); oR = eas_func_dwn(xR, fade, pfad);
                    pfad--; if(pfad<0)spprev=sp;
                }else{ oL = oR = 0.; spprev=sp; }  //.........<-Everything Off
            }
        }
        else//LOOPER_(not_oneshot)_SECTION
        {       //................................LOOPER-RECORDING
            if((rs!=rspprev)&&(rfad<0))
            {
                if(((rspprev==0)&&(rs!=0))&&(rfad<0)){ rfad=fade-1; wram=1; rspprev=rs; }
                else if (((rspprev!=0)&&(rs==0))&&(rfad<0)){ rfad=fade-1; }
                else if(((rspprev<0)&&(rs>0))||((rspprev>0)&&(rs<0))){ rxfad=fade; rxpos=rpos; }
            }
            
            if(rs!=0)
            {
                if((rinpos!=rphprev)&&(rxfad<0))
                {                                 //crossfade positional changes
                    if(rwr){ if((rinpos>nrend)&&(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    else{ if((rinpos>nrend)||(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    if(rfad<0){rxfad=fade; rxpos=rpos;} rphprev=rpos=rinpos;
                }
                
                if(rparamch)
                {                                 //window changes
                    if(rwr){
                        if(((rpos>rend)&&(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }else{
                        if(((rpos>rend)||(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }
                }
                                            //recording-crossfade at loop boundaries
                if((rinOderAus)&&(rxfad<0))
                {                                //rparamch = flag to notify when there's incoming..
                    if(rparamch){ rparamch=0; nrend=rend; nrstart=rstart; }//..message-set loop-bounds
                    if(rwr){
                        if((rpos>nrend)&&(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }else{
                        if((rpos>nrend)||(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }
                }
                else        //if starting off outside of loop-boundaries, record until within...
                {       //..and once within, save messages(rstart/rend) internally(nrstart/nrend)
                    if(rwr)
                    {if((rpos<=rend)||(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                    else
                    {if((rpos<=rend)&&(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rxfad>=0) //....Crossfades(happen at loop-points and phase changes)
                {                    //crossfade-point buffer~-boundary constraint and increment
                    if(rxpos>bframes)rxpos=0; if(rxpos<0)rxpos=bframes; rxdx=trunc(rxpos); rxpos+=rs;
                    b[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc],fade,1-odb,rxfad);
                    b[rxdx*bnc+1]=eas_func_dwn(rcR,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc+1],fade,1-odb,rxfad);
                    b[rndx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rxfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rxfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rxfad);
                    rxfad--;
                }
                else if(rfad>=0)     //...................Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rfad);rfad--;
                } else { b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); b[rndx*bnc+1]=rcR+(b[rndx*bnc+1]*odb); }
            }
            else
            {       //...................................Fade-Out
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_dwn(rcR,fade,rfad)+eas_rec_dwn(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; rspprev=rs; }
                }else{ rspprev=rs; }
            }
            
                //.....................................LOOPER-PLAYBACK
            if((sp!=spprev)&&(xfad<0))
            {if((((spprev==0.0)&&(sp!=0.0))||((spprev!=0.0)&&(sp==0.0)))&&(pfad<0)){pfad=fade;}//onoff
            else if(((spprev<0)&&(sp>0))||((spprev>0)&&(sp<0))){xfad=fade;xpos=dpos;}}//xfad neg<->pos
            
            if(sp!=0.0)
            {
                if((inpos!=phprev)&&(xfad<0))
                {                                       //crossfade positional changes
                    if(wr){ if((inpos>nend)&&(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    else{ if((inpos>nend)||(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    if(pfad<0){xfad=fade; xpos=dpos;} phprev=dpos=inpos;
                }
                if(pparamch)
                {                                      //window changes
                    if(wr){
                        if(((dpos>pend)&&(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }else{
                        if(((dpos>pend)||(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }
                }
                                                //crossfade at playback-loop boundaries
                if((inOderAus)&&(xfad<0))
                {
                    if(pparamch){ pparamch=0; nend=pend; nstart=pstart; }
                    if(wr){
                        if((dpos>nend)&&(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }else{
                        if((dpos>nend)||(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }
                }
                else        //if starting off outside of loop-boundaries, play until within...
                {         //..and once within, save messages(pstart/pend) internally(nstart/nend)
                    if(wr)
                    {if((dpos<=pend)||(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                    else
                    {if((dpos<=pend)&&(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);

                if(pfad>=0)   //.......................Fade-In
                { oL=eas_func_up(xL,fade,pfad); oR=eas_func_up(xR,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);
                    oR = eas_func_up(xR, fade, xfad) + eas_func_dwn(oR, fade, xfad); xfad--;
                }else{ oL=xL; oR=xR; }  //..<-Regular Playback
                spprev=sp;
            }
            else
            {  //......................................Fade-Out
                if(pfad>=0)
                {
                    if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                    if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                    dpos+=spprev;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_dwn(xL, fade, pfad); oR = eas_func_dwn(xR, fade, pfad); pfad--;
                    if(pfad<0)spprev=sp;
                }else{ oL = oR = 0.; spprev=sp; }  //<-Everything Off
            }
            
            if(wram)    //..calculate abs-difference between playhead and recordhead..
            {   //..if difference is within fade time, use difference to drive ducking of playhead
                if(sp!=(double)rs)
                {
                    dfrnc = fabs(rpos-dpos);
                    if(dfrnc<(fade*2.7489))
                    { oL=eas_func_dwn(oL,(fade*2.7489),dfrnc); oR=eas_func_dwn(oR,(fade*2.7489),dfrnc); }
                }
            }
        }
        
        oPh=dpos/bframes; oRPh=rpos/bframes;//<-convert samp-index 2 phase
        *outL++ = oL; *outR++ = oR; *outPh++ = oPh; *outRPh++ = oRPh;
    }
    if(wram)buffer_setdirty(buffer); buffer_unlocksamples(buffer); x->nshot=nshot; x->rnshot=rnshot;
    x->pparamchange=pparamch; x->rparamchange=rparamch; x->xfad=xfad; x->pfad=pfad; x->xpos=xpos;
    x->phprev=phprev; x->rphprev=rphprev; x->spprev=spprev; x->rspprev=rspprev; x->rxpos=rxpos;
    x->dpos=dpos; x->rpos=rpos; x->inOderAus=inOderAus; x->rinOderAus=rinOderAus; x->rfad=rfad;
    x->nend=nend; x->nstart=nstart; x->nrend=nrend; x->nrstart=nrstart; x->rxfad=rxfad; x->wram=wram;
    return;
    
zero: while(n--){ *outL++=0.; *outR++=0.; *outPh++=0.; *outRPh++=0.; }
}

void rez_psperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
    t_double        *ph = is[0];
    t_double        *spd = is[1];
    t_double        *outL = os[0];
    t_double        *outR = os[1];
    t_double        *outPh = os[2];

    long         n=smps;                                           float *b;
    t_double     phprev,inpos,dpos,sp,frac;
    t_double     spprev,oL,oR,oPh,xL,xR,xpos,scldsr;
    t_ptr_int    bnc,bframes,pstart,pend,nstart,nend,indx,indxz,indxb,indxc,xfad,pfad,fade;
    t_bool       wr,pparamch,inOderAus,oneshot;
    t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

    b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot; nstart=x->nstart;
    bnc=x->bchan; bframes=x->bframes-1; pstart=x->pstart; pend=x->pend; inOderAus=x->inOderAus;
    pparamch=x->pparamchange; pfad=x->pfad; xfad=x->xfad; fade=x->fade; xpos=x->xpos;
    wr=x->wr; phprev=x->phprev; spprev=x->spprev; dpos=x->dpos; scldsr=x->srscale; nend=x->nend;
    
    while(n--)
    {
        inpos=*ph++; sp=*spd++;
        inpos=(inpos<0.0)?0.0:((inpos>1.0)?1.0:inpos); inpos*=bframes; sp*=scldsr;

        if(oneshot)
        {
            if((spprev==0.0)&&(sp!=0.0))
            { dpos=inpos; spprev=sp; pfad=fade; inOderAus=1; nend=pend; nstart=pstart; }
            else if(((spprev!=0.0)&&(sp==0.0))&&(inOderAus==1)){ pfad=fade; inOderAus=0; }
            if (inOderAus)
            {
                if(wr){
                    if(((dpos>(nend-(fade+1)))&&(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }else
                {
                    if(((dpos>(nend-(fade+1)))||(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);

                if(pfad>=0)   //.......................Fade-In
                { oL=eas_func_up(xL,fade,pfad); oR=eas_func_up(xR,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);
                    oR = eas_func_up(xR, fade, xfad) + eas_func_dwn(oR, fade, xfad); xfad--;
                }else{ oL=xL; oR=xR; }  //..<-Regular Playback + Phase-History
            }
            else
            {//......................................Fade-Out
                    if(pfad>=0)
                    {
                        if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                        if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                        dpos+=spprev;
                        interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                        xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                        xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                        oL = eas_func_dwn(xL, fade, pfad); oR = eas_func_dwn(xR, fade, pfad); pfad--;
                    }else{ oL = oR = 0.; spprev=sp; }  //.........<-Everything Off
            }
        }
        else
        {
            if((sp!=spprev)&&(xfad<0))
            {if((((spprev==0.0)&&(sp!=0.0))||((spprev!=0.0)&&(sp==0.0)))&&(pfad<0)){pfad=fade;}//onoff
            else if(((spprev<0)&&(sp>0))||((spprev>0)&&(sp<0))){xfad=fade;xpos=dpos;}}//xfad neg<->pos
            
            if(sp!=0.0)
            {
                if((inpos!=phprev)&&(xfad<0))
                {                                       //crossfade positional changes
                    if(wr){ if((inpos>nend)&&(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    else{ if((inpos>nend)||(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    if(pfad<0){xfad=fade; xpos=dpos;} phprev=dpos=inpos;
                }
                if(pparamch)
                {                                               //window changes
                    if(wr){
                        if(((dpos>pend)&&(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }else{
                        if(((dpos>pend)||(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }
                }
                
                if((inOderAus)&&(xfad<0))
                {
                    if(pparamch){ pparamch=0; nend=pend; nstart=pstart; }
                    if(wr){
                        if((dpos>nend)&&(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }else{
                        if((dpos>nend)||(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }
                }
                else
                {
                    if(wr)
                    {if((dpos<=pend)||(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                    else
                    {if((dpos<=pend)&&(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);

                if(pfad>=0)   //.......................Fade-In
                { oL=eas_func_up(xL,fade,pfad); oR=eas_func_up(xR,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);
                    oR = eas_func_up(xR, fade, xfad) + eas_func_dwn(oR, fade, xfad); xfad--;
                }else{ oL=xL; oR=xR; }  //..<-Regular Playback + Phase-History
                spprev=sp;
            }
            else
            {  //......................................Fade-Out
                if(pfad>=0)
                {
                    if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                    if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                    dpos+=spprev;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    xR = interp(frac, b[indxz*bnc+1], b[indx*bnc+1], b[indxb*bnc+1], b[indxc*bnc+1]);
                    oL = eas_func_dwn(xL, fade, pfad); oR = eas_func_dwn(xR, fade, pfad); pfad--;
                    if(pfad<=0)spprev=sp;
                }else{ oL = oR = 0.; }  //.........<-Everything Off
            }
        }
        
        oPh=dpos/bframes; //<-convert samp-index 2 phase
        *outL++ = oL; *outR++ = oR; *outPh++ = oPh;
    }
    buffer_unlocksamples(buffer);
    x->pparamchange=pparamch; x->xfad=xfad; x->pfad=pfad; x->xpos=xpos; x->phprev=phprev;
    x->dpos=dpos; x->inOderAus=inOderAus; x->nend=nend; x->nstart=nstart; x->spprev=spprev;
    return;
    
zero: while(n--){ *outL++=0.; *outR++=0.; *outPh++=0.; }
}

void rez_rsperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
    t_double        *rph = is[0];
    t_double        *rspd = is[1];
    t_double        *ovdb = is[2];
    t_double        *rinL = is[3];
    t_double        *rinR = is[4];
    t_double        *outRPh = os[0];

    long         n=smps;                                           float *b;
    t_double     rcL,rcR,rpos,rxpos,rsp,oRPh,odb,rphprev,rinpos;
    t_ptr_int    rndx,rxdx,bnc,bframes,rstart,rend,nrstart,nrend,rfad,rxfad,fade,rs,rspprev;
    t_bool       rwr,rparamch,rinOderAus,wram,oneshot;
    t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

    b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot; wram=x->wram;
    bnc=x->bchan; bframes=x->bframes-1; fade=x->fade; nrend=x->nrend; nrstart=x->nrstart; rwr=x->rwr;
    rparamch=x->rparamchange; rfad=x->rfad; rstart=x->rstart; rend=x->rend; rphprev=x->rphprev;
    rxfad=x->rxfad; rxpos=x->rxpos; rspprev=x->rspprev; rinOderAus=x->rinOderAus; rpos=x->rpos;
    
    while(n--)
    {
        rinpos=*rph++; rsp=*rspd++; odb=*ovdb++; rcL=*rinL++; rcR=*rinR++;
        rinpos=(rinpos<0.0)?0.0:((rinpos>1.0)?1.0:rinpos); rinpos*=bframes;
        if(rsp>0.0) rs=1; else if(rsp<0.0) rs=-1; else rs=0;
        
                        //.........................................RECORDING
        if(oneshot)
        {
            if((rspprev==0)&&(rs!=0))
            { rpos=rinpos; rspprev=rs; rfad=fade; wram=1; rinOderAus=1; nrend=rend; nrstart=rstart; }
            else if(((rspprev!=0)&&(rs==0))&&(rinOderAus==1)){ rfad=fade; rinOderAus=0; }
            if (rinOderAus)
            {
                if(rwr){
                    if(((rpos>(nrend-(fade+1)))&&(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }else
                {
                    if(((rpos>(nrend-(fade+1)))||(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rfad>=0)     //........................Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--;
                }else{b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); b[rndx*bnc+1]=rcR+(b[rndx*bnc+1]*odb);};
            }
            else
            {
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_dwn(rcR,fade,rfad)+eas_rec_dwn(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; }
                }else rspprev=rs;
            }
        }
        else
        {
            if((rs!=rspprev)&&(rfad<0))
            {
                if(((rspprev==0)&&(rs!=0))&&(rfad<0)){ rfad=fade; wram=1; rspprev=rs; }
                else if (((rspprev!=0)&&(rs==0))&&(rfad<0)){ rfad=fade; }
                else if(((rspprev<0)&&(rs>0))||((rspprev>0)&&(rs<0))){rxfad=fade;rxpos=rpos;}
            }//onoff
            
            if(rs!=0)
            {
                if((rinpos!=rphprev)&&(rxfad<0))
                {                                 //crossfade positional changes
                    if(rwr){ if((rinpos>nrend)&&(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    else{ if((rinpos>nrend)||(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    if(rfad<0){rxfad=fade; rxpos=rpos;} rphprev=rpos=rinpos;
                }
                
                if(rparamch)
                {                                 //window changes
                    if(rwr){
                        if(((rpos>rend)&&(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }else{
                        if(((rpos>rend)||(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }
                }
                
                if((rinOderAus)&&(rxfad<0))
                {
                    if(rparamch){ rparamch=0; nrend=rend; nrstart=rstart; }
                    if(rwr){
                        if((rpos>nrend)&&(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }else{
                        if((rpos>nrend)||(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }
                }
                else
                {
                    if(rwr)
                    {if((rpos<=rend)||(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                    else
                    {if((rpos<=rend)&&(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rxfad>=0) //....Crossfades(happen at loop-points and phase changes)
                {
                    if(rxpos>bframes)rxpos=0; if(rxpos<0)rxpos=bframes; rxdx=trunc(rxpos); rxpos+=rs;
                    b[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc],fade,1-odb,rxfad);
                    b[rxdx*bnc+1]=eas_func_dwn(rcR,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc+1],fade,1-odb,rxfad);
                    b[rndx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rxfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rxfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rxfad);
                    rxfad--;
                }
                else if(rfad>=0)     //........................Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_up(rcR,fade,rfad)+eas_rec_up(b[rndx*bnc+1],fade,1-odb,rfad);rfad--;
                } else { b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); b[rndx*bnc+1]=rcR+(b[rndx*bnc+1]*odb); }
            }
            else
            {       //.......................................Fade-Out
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    b[rndx*bnc+1]=eas_func_dwn(rcR,fade,rfad)+eas_rec_dwn(b[rndx*bnc+1],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; rspprev=rs; }
                }
            }
        }
        oRPh=rpos/bframes;//<-convert samp-index 2 phase
        *outRPh++ = oRPh;
    }
    if(wram)buffer_setdirty(buffer); buffer_unlocksamples(buffer); x->rparamchange=rparamch;
    x->rphprev=rphprev; x->rspprev=rspprev; x->rxpos=rxpos; x->rpos=rpos; x->rinOderAus=rinOderAus;
    x->rfad=rfad; x->nrend=nrend; x->nrstart=nrstart; x->rxfad=rxfad; x->wram=wram;
    return;
    
zero: while(n--){ *outRPh++=0.; }
}

void rez_mperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
    t_double        *ph = is[0];
    t_double        *spd = is[1];
    t_double        *rph = is[2];
    t_double        *rspd = is[3];
    t_double        *ovdb = is[4];
    t_double        *rinL = is[5];
    t_double        *outL = os[0];
    t_double        *outPh = os[1];
    t_double        *outRPh = os[2];

    long         n=smps;                                           float *b;
    t_double     phprev,inpos,dpos,sp,rcL,frac,dfrnc,rpos,rxpos,rsp;
    t_double     spprev,oL,oPh,oRPh,xL,xpos,scldsr,odb,rphprev,rinpos;
    t_ptr_int    rndx,rxdx,bnc,bframes,pstart,pend,rstart,rend,nstart,nend,nrstart,nrend;
    t_ptr_int    indx,indxz,indxb,indxc,xfad,pfad,rfad,rxfad,fade,rs,rspprev;
    t_bool       wr,rwr,pparamch,rparamch,inOderAus,rinOderAus,wram,oneshot;
    t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

    b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot;
    bnc=x->bchan; bframes=x->bframes-1; pstart=x->pstart; pend=x->pend; inOderAus=x->inOderAus;
    pparamch=x->pparamchange; pfad=x->pfad; xfad=x->xfad; fade=x->fade; xpos=x->xpos; rpos=x->rpos;
    wr=x->wr;  rparamch=x->rparamchange; rfad=x->rfad; rstart=x->rstart; rend=x->rend;
    phprev=x->phprev; rphprev=x->rphprev; spprev=x->spprev;  dpos=x->dpos; wram=x->wram;
    scldsr=x->srscale; nend=x->nend; nstart=x->nstart; nrend=x->nrend; nrstart=x->nrstart;
    rxfad=x->rxfad; rxpos=x->rxpos; rspprev=x->rspprev; rwr=x->rwr; rinOderAus=x->rinOderAus;
    
    while(n--)
    {
        inpos=*ph++; sp=*spd++; rinpos=*rph++; rsp=*rspd++; odb=*ovdb++; rcL=*rinL++;
        inpos=(inpos<0.0)?0.0:((inpos>1.0)?1.0:inpos); inpos*=bframes; sp*=scldsr;
        rinpos=(rinpos<0.0)?0.0:((rinpos>1.0)?1.0:rinpos); rinpos*=bframes;
        if(rsp>0.0) rs=1; else if(rsp<0.0) rs=-1; else rs=0;
        
                        //.........................................RECORDING
        if(oneshot)
        {
            if((rspprev==0)&&(rs!=0))
            { rpos=rinpos; rspprev=rs; rfad=fade; wram=1; rinOderAus=1; nrend=rend; nrstart=rstart; }
            else if(((rspprev!=0)&&(rs==0))&&(rinOderAus==1)){ rfad=fade; rinOderAus=0; }
            if (rinOderAus)
            {
                if(rwr){
                    if(((rpos>(nrend-(fade+1)))&&(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }else
                {
                    if(((rpos>(nrend-(fade+1)))||(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rfad>=0)     //........................Fade-In
                { b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);rfad--; }
                else{ b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); };
            }
            else
            {
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; }
                }else rspprev=rs;
            }
            
            //.........................................PLAYBACK
            if((spprev==0.0)&&(sp!=0.0))
            { dpos=inpos; spprev=sp; pfad=fade; inOderAus=1; nend=pend; nstart=pstart; }
            else if(((spprev!=0.0)&&(sp==0.0))&&(inOderAus==1)){ pfad=fade; inOderAus=0; }
            if (inOderAus)
            {
                if(wr){
                    if(((dpos>(nend-(fade+1)))&&(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }else
                {
                    if(((dpos>(nend-(fade+1)))||(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                //.......................Fade-In
                if(pfad>=0){ oL=eas_func_up(xL,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad); xfad--;
                }else{ oL=xL; }  //..<-Regular Playback + Phase-History
            }
            else
            {//......................................Fade-Out
                    if(pfad>=0)
                    {
                        if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                        if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                        dpos+=spprev;
                        interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                        xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                        oL = eas_func_dwn(xL, fade, pfad); pfad--;
                    }else{ oL = 0.; spprev=sp; }  //.........<-Everything Off
            }
        }
        else
        {
            if((rs!=rspprev)&&(rfad<0))
            {
                if(((rspprev==0)&&(rs!=0))&&(rfad<0)){ rfad=fade; wram=1; rspprev=rs; }
                else if (((rspprev!=0)&&(rs==0))&&(rfad<0)){ rfad=fade; }
                else if(((rspprev<0)&&(rs>0))||((rspprev>0)&&(rs<0))){rxfad=fade;rxpos=rpos;}
            }//onoff
            
            if(rs!=0)
            {
                if((rinpos!=rphprev)&&(rxfad<0))
                {                                 //crossfade positional changes
                    if(rwr){ if((rinpos>nrend)&&(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    else{ if((rinpos>nrend)||(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    if(rfad<0){rxfad=fade; rxpos=rpos;} rphprev=rpos=rinpos;
                }
                
                if(rparamch)
                {                                 //window changes
                    if(rwr){
                        if(((rpos>rend)&&(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }else{
                        if(((rpos>rend)||(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }
                }
                
                if((rinOderAus)&&(rxfad<0))
                {
                    if(rparamch){ rparamch=0; nrend=rend; nrstart=rstart; }
                    if(rwr){
                        if((rpos>nrend)&&(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }else{
                        if((rpos>nrend)||(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }
                }
                else
                {
                    if(rwr)
                    {if((rpos<=rend)||(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                    else
                    {if((rpos<=rend)&&(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rxfad>=0) //....Crossfades(happen at loop-points and phase changes)
                {
                    if(rxpos>bframes)rxpos=0; if(rxpos<0)rxpos=bframes; rxdx=trunc(rxpos); rxpos+=rs;
                    b[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc],fade,1-odb,rxfad);
                    b[rndx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rxfad);
                    rxfad--;
                }
                else if(rfad>=0)     //........................Fade-In
                { b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad); }
                else{ b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); }
            }
            else
            {       //.......................................Fade-Out
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; rspprev=rs; }
                }
            }
            
            //.........................................PLAYBACK
            if((sp!=spprev)&&(xfad<0))
            {if((((spprev==0.0)&&(sp!=0.0))||((spprev!=0.0)&&(sp==0.0)))&&(pfad<0)){pfad=fade;}//onoff
            else if(((spprev<0)&&(sp>0))||((spprev>0)&&(sp<0))){xfad=fade;xpos=dpos;}}//xfad neg<->pos
            
            if(sp!=0.0)
            {
                if((inpos!=phprev)&&(xfad<0))
                {                                       //crossfade positional changes
                    if(wr){ if((inpos>nend)&&(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    else{ if((inpos>nend)||(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    if(pfad<0){xfad=fade; xpos=dpos;} phprev=dpos=inpos;
                }
                if(pparamch)
                {                                               //window changes
                    if(wr){
                        if(((dpos>pend)&&(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }else{
                        if(((dpos>pend)||(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }
                }
                
                if((inOderAus)&&(xfad<0))
                {
                    if(pparamch){ pparamch=0; nend=pend; nstart=pstart; }
                    if(wr){
                        if((dpos>nend)&&(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }else{
                        if((dpos>nend)||(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }
                }
                else
                {
                    if(wr)
                    {if((dpos<=pend)||(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                    else
                    {if((dpos<=pend)&&(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                //.......................Fade-In
                if(pfad>=0){ oL=eas_func_up(xL,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad); xfad--;
                }else{ oL=xL; }  //..<-Regular Playback + Phase-History
                spprev=sp;
            }
            else
            {  //......................................Fade-Out
                if(pfad>=0)
                {
                    if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                    if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                    dpos+=spprev;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_dwn(xL, fade, pfad); pfad--;
                    if(pfad<=0)spprev=sp;
                }else{ oL = 0.; }  //.........<-Everything Off
            }
            
            if(wram)
            {
                if(sp!=(double)rs)
                { dfrnc=fabs(rpos-dpos); if(dfrnc<=(fade*2)){ oL=eas_func_dwn(oL,(fade*2),dfrnc); } }
            }
        }
        
        oPh=dpos/bframes; oRPh=rpos/bframes;//<-convert samp-index 2 phase
        *outL++ = oL; *outPh++ = oPh; *outRPh++ = oRPh;
    }
    if(wram)buffer_setdirty(buffer); buffer_unlocksamples(buffer);
    x->pparamchange=pparamch; x->rparamchange=rparamch; x->xfad=xfad; x->pfad=pfad; x->xpos=xpos;
    x->phprev=phprev; x->rphprev=rphprev; x->spprev=spprev; x->rspprev=rspprev; x->rxpos=rxpos;
    x->dpos=dpos; x->rpos=rpos; x->inOderAus=inOderAus; x->rinOderAus=rinOderAus; x->rfad=rfad;
    x->nend=nend; x->nstart=nstart; x->nrend=nrend; x->nrstart=nrstart; x->rxfad=rxfad; x->wram=wram;
    return;
    
zero: while(n--){ *outL++=0.; *outPh++=0.; *outRPh++=0.; }
}

void rez_pmperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
    t_double        *ph = is[0];
    t_double        *spd = is[1];
    t_double        *outL = os[0];
    t_double        *outPh = os[1];

    long         n=smps;                                           float *b;
    t_double     phprev,inpos,dpos,sp,frac;
    t_double     spprev,oL,oPh,xL,xpos,scldsr;
    t_ptr_int    bnc,bframes,pstart,pend,nstart,nend,indx,indxz,indxb,indxc,xfad,pfad,fade;
    t_bool       wr,pparamch,inOderAus,oneshot;
    t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

    b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot; nstart=x->nstart;
    bnc=x->bchan; bframes=x->bframes-1; pstart=x->pstart; pend=x->pend; inOderAus=x->inOderAus;
    pparamch=x->pparamchange; pfad=x->pfad; xfad=x->xfad; fade=x->fade; xpos=x->xpos;
    wr=x->wr; phprev=x->phprev; spprev=x->spprev; dpos=x->dpos; scldsr=x->srscale; nend=x->nend;
    
    while(n--)
    {
        inpos=*ph++; sp=*spd++;
        inpos=(inpos<0.0)?0.0:((inpos>1.0)?1.0:inpos); inpos*=bframes; sp*=scldsr;

        if(oneshot)
        {
            if((spprev==0.0)&&(sp!=0.0))
            { dpos=inpos; spprev=sp; pfad=fade; inOderAus=1; nend=pend; nstart=pstart; }
            else if(((spprev!=0.0)&&(sp==0.0))&&(inOderAus==1)){ pfad=fade; inOderAus=0; }
            if (inOderAus)
            {
                if(wr){
                    if(((dpos>(nend-(fade+1)))&&(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }else
                {
                    if(((dpos>(nend-(fade+1)))||(dpos<(nstart+(fade+1))))&&(pfad<0))
                    { inOderAus=0; pfad=fade; }
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                //.......................Fade-In
                if(pfad>=0){ oL=eas_func_up(xL,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad);
                }else{ oL=xL; }  //..<-Regular Playback + Phase-History
            }
            else
            {//......................................Fade-Out
                    if(pfad>=0)
                    {
                        if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                        if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                        dpos+=spprev;
                        interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                        xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                        oL = eas_func_dwn(xL, fade, pfad); pfad--;
                    }else{ oL = 0.; spprev=sp; }  //.........<-Everything Off
            }
        }
        else
        {
            if((sp!=spprev)&&(xfad<0))
            {if((((spprev==0.0)&&(sp!=0.0))||((spprev!=0.0)&&(sp==0.0)))&&(pfad<0)){pfad=fade;}//onoff
            else if(((spprev<0)&&(sp>0))||((spprev>0)&&(sp<0))){xfad=fade;xpos=dpos;}}//xfad neg<->pos
            
            if(sp!=0.0)
            {
                if((inpos!=phprev)&&(xfad<0))
                {                                       //crossfade positional changes
                    if(wr){ if((inpos>nend)&&(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    else{ if((inpos>nend)||(inpos<nstart))inOderAus=0; else inOderAus=1; }
                    if(pfad<0){xfad=fade; xpos=dpos;} phprev=dpos=inpos;
                }
                if(pparamch)
                {                                               //window changes
                    if(wr){
                        if(((dpos>pend)&&(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }else{
                        if(((dpos>pend)||(dpos<pstart))&&((pend>nend)||(pstart<nstart)))
                            inOderAus=0; else inOderAus=1;
                    }
                }
                
                if((inOderAus)&&(xfad<0))
                {
                    if(pparamch){ pparamch=0; nend=pend; nstart=pstart; }
                    if(wr){
                        if((dpos>nend)&&(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }else{
                        if((dpos>nend)||(dpos<nstart))
                        {
                            xfad=fade; xpos=dpos;
                            if(sp>0.0) dpos=(dpos-nend)+nstart; else dpos=nend-(nstart-dpos);
                        }
                    }
                }
                else
                {
                    if(wr)
                    {if((dpos<=pend)||(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                    else
                    {if((dpos<=pend)&&(dpos>=pstart)){pparamch=0; inOderAus=1; nend=pend; nstart=pstart;}}
                }
                
                if(dpos>bframes)dpos-=(bframes+1);  if(dpos<0)dpos+=(bframes+1);  indx=trunc(dpos);
                if(sp>0){ frac=dpos-indx; }else if(sp<0){ frac=1.0-(dpos-indx); }else frac=0.0; dpos+=sp;
                interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);

                if(pfad>=0)   //.......................Fade-In
                { oL=eas_func_up(xL,fade,pfad); pfad--; }
                else if(xfad>=0)//..Crossfades(happen at loop-points,phase-changes,and abrupt pos/neg)
                {
                    if(xpos>bframes)xpos-=(bframes+1); if(xpos<0)xpos+=(bframes+1); indx=trunc(xpos);
                    if(sp>0){frac=xpos-indx;}else if(sp<0){frac=1.0-(xpos-indx);}else frac=0.0;
                    xpos+=sp;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    oL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_up(xL, fade, xfad) + eas_func_dwn(oL, fade, xfad); xfad--;
                }else{ oL=xL; }  //..<-Regular Playback + Phase-History
                spprev=sp;
            }
            else
            {  //......................................Fade-Out
                if(pfad>=0)
                {
                    if(dpos>bframes)dpos-=(bframes+1); if(dpos<0)dpos+=(bframes+1); indx=trunc(dpos);
                    if(spd>0){frac=dpos-indx;}else if(spd<0){frac=1.0-(dpos-indx);}else frac=0.0;
                    dpos+=spprev;
                    interp_index(indx,&indxz,&indxb,&indxc,sp,bframes);
                    xL = interp(frac, b[indxz*bnc], b[indx*bnc], b[indxb*bnc], b[indxc*bnc]);
                    oL = eas_func_dwn(xL, fade, pfad); pfad--;
                    if(pfad<=0)spprev=sp;
                }else{ oL = 0.; }  //.........<-Everything Off
            }
        }
        
        oPh=dpos/bframes; //<-convert samp-index 2 phase
        *outL++ = oL; *outPh++ = oPh;
    }
    buffer_unlocksamples(buffer);
    x->pparamchange=pparamch; x->xfad=xfad; x->pfad=pfad; x->xpos=xpos; x->phprev=phprev;
    x->dpos=dpos; x->inOderAus=inOderAus; x->nend=nend; x->nstart=nstart; x->spprev=spprev;
    return;
    
zero: while(n--){ *outL++=0.; *outPh++=0.; }
}

void rez_rmperform64(t_rez *x, t_object *d64, double **is, long numis, double **os, long numos, long smps, long flgs, void *usrp)
{
    t_double        *rph = is[0];
    t_double        *rspd = is[1];
    t_double        *ovdb = is[2];
    t_double        *rinL = is[3];
    t_double        *outRPh = os[0];

    long         n=smps;                                           float *b;
    t_double     rcL,rpos,rxpos,rsp,oRPh,odb,rphprev,rinpos;
    t_ptr_int    rndx,rxdx,bnc,bframes,rstart,rend,nrstart,nrend,rfad,rxfad,fade,rs,rspprev;
    t_bool       rwr,rparamch,rinOderAus,wram,oneshot;
    t_buffer_obj *buffer=buffer_ref_getobject(x->buf);

    b = buffer_locksamples(buffer);                   if(!b||x->obj.z_disabled) goto zero;
    if(x->buf_modd){ x->buf_modd=false; rez_limits(x); } oneshot=x->oneshot; wram=x->wram;
    bnc=x->bchan; bframes=x->bframes-1; fade=x->fade; nrend=x->nrend; nrstart=x->nrstart; rwr=x->rwr;
    rparamch=x->rparamchange; rfad=x->rfad; rstart=x->rstart; rend=x->rend; rphprev=x->rphprev;
    rxfad=x->rxfad; rxpos=x->rxpos; rspprev=x->rspprev; rinOderAus=x->rinOderAus; rpos=x->rpos;
    
    while(n--)
    {
        rinpos=*rph++; rsp=*rspd++; odb=*ovdb++; rcL=*rinL++;
        rinpos=(rinpos<0.0)?0.0:((rinpos>1.0)?1.0:rinpos); rinpos*=bframes;
        if(rsp>0.0) rs=1; else if(rsp<0.0) rs=-1; else rs=0;
        
                        //.........................................RECORDING
        if(oneshot)
        {
            if((rspprev==0)&&(rs!=0))
            { rpos=rinpos; rspprev=rs; rfad=fade; wram=1; rinOderAus=1; nrend=rend; nrstart=rstart; }
            else if(((rspprev!=0)&&(rs==0))&&(rinOderAus==1)){ rfad=fade; rinOderAus=0; }
            if (rinOderAus)
            {
                if(rwr){
                    if(((rpos>(nrend-(fade+1)))&&(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }else
                {
                    if(((rpos>(nrend-(fade+1)))||(rpos<(nrstart+(fade+1))))&&(rfad<0))
                    { rinOderAus=0; rfad=fade; }
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rfad>=0)     //........................Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--;
                }else{ b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); };
            }
            else
            {
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; }
                }else rspprev=rs;
            }
        }
        else
        {
            if((rs!=rspprev)&&(rfad<0))
            {
                if(((rspprev==0)&&(rs!=0))&&(rfad<0)){ rfad=fade; wram=1; rspprev=rs; }
                else if (((rspprev!=0)&&(rs==0))&&(rfad<0)){ rfad=fade; }
                else if(((rspprev<0)&&(rs>0))||((rspprev>0)&&(rs<0))){rxfad=fade;rxpos=rpos;}
            }//onoff
            
            if(rs!=0)
            {
                if((rinpos!=rphprev)&&(rxfad<0))
                {                                 //crossfade positional changes
                    if(rwr){ if((rinpos>nrend)&&(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    else{ if((rinpos>nrend)||(rinpos<nrstart))rinOderAus=0; else rinOderAus=1; }
                    if(rfad<0){rxfad=fade; rxpos=rpos;} rphprev=rpos=rinpos;
                }
                
                if(rparamch)
                {                                 //window changes
                    if(rwr){
                        if(((rpos>rend)&&(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }else{
                        if(((rpos>rend)||(rpos<rstart))&&((rend>nrend)||(rstart<nrstart)))
                            rinOderAus=0; else rinOderAus=1;
                    }
                }
                
                if((rinOderAus)&&(rxfad<0))
                {
                    if(rparamch){ rparamch=0; nrend=rend; nrstart=rstart; }
                    if(rwr){
                        if((rpos>nrend)&&(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }else{
                        if((rpos>nrend)||(rpos<nrstart))
                        {
                            rxfad=fade;
                            if(rs>0){ rpos=nrstart; rxpos=nrend; }else{ rpos=nrend; rxpos=nrstart;}
                        }
                    }
                }
                else
                {
                    if(rwr)
                    {if((rpos<=rend)||(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                    else
                    {if((rpos<=rend)&&(rpos>=rstart))
                    {rparamch=0; rinOderAus=1; nrend=rend; nrstart=rstart;}}
                }
                
                if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rs;
                if(rxfad>=0) //....Crossfades(happen at loop-points and phase changes)
                {
                    if(rxpos>bframes)rxpos=0; if(rxpos<0)rxpos=bframes; rxdx=trunc(rxpos); rxpos+=rs;
                    b[rxdx*bnc]=eas_func_dwn(rcL,fade,rxfad)+eas_rec_dwn(b[rxdx*bnc],fade,1-odb,rxfad);
                    b[rndx*bnc]=eas_func_up(rcL,fade,rxfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rxfad);
                    rxfad--;
                }
                else if(rfad>=0)     //........................Fade-In
                {
                    b[rndx*bnc]=eas_func_up(rcL,fade,rfad)+eas_rec_up(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--;
                } else { b[rndx*bnc]=rcL+(b[rndx*bnc]*odb); }
            }
            else
            {       //.......................................Fade-Out
                if(rfad>=0)
                {
                    if(rpos>bframes)rpos=0;  if(rpos<0)rpos=bframes; rndx=trunc(rpos); rpos+=rspprev;
                    b[rndx*bnc]=eas_func_dwn(rcL,fade,rfad)+eas_rec_dwn(b[rndx*bnc],fade,1-odb,rfad);
                    rfad--; if(rfad<0) { wram=0; rspprev=rs; }
                }
            }
        }
        oRPh=rpos/bframes;//<-convert samp-index 2 phase
        *outRPh++ = oRPh;
    }
    if(wram)buffer_setdirty(buffer); buffer_unlocksamples(buffer); x->rparamchange=rparamch;
    x->rphprev=rphprev; x->rspprev=rspprev; x->rxpos=rxpos; x->rpos=rpos; x->rinOderAus=rinOderAus;
    x->rfad=rfad; x->nrend=nrend; x->nrstart=nrstart; x->rxfad=rxfad; x->wram=wram;
    return;
    
zero: while(n--){ *outRPh++=0.; }
}
