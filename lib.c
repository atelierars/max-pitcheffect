#include"ext.h"            // standard Max include, always required (except in Jitter)
#include"ext_obex.h"        // required for "new" style objects
#include"z_dsp.h"            // required for MSP objects
#include<Accelerate/Accelerate.h>
#include<simd/simd.h>
typedef struct {
    t_pxobject super;
	double const pitch;
	double const width;
	double const verge;
	double const level;
	intptr_t frame;
	intptr_t index;
	intptr_t count;
	intptr_t phase;
	double * cache;
	vDSP_DFT_SetupD forward;
	vDSP_DFT_SetupD inverse;
} t_pitcheffect;
C74_HIDDEN void hilbert(double * const target, double const lambda, uintptr_t const length) {
	int const center = (int const) length / 2 + 1;
	vDSP_vrampD((double const[]){0}, (double const[]){1}, target, 1, center);
	vDSP_vsqD(target, 1, target, 1, center);
	vDSP_vsdivD(target, 1, (double const[]){center*center/-(isnormal(lambda)?lambda*lambda:1)}, target, 1, center);
	vvexp(target, target, &center);
	vDSP_vsmulD(target + 1, 1, (double const[]){2}, target + 1, 1, (length - 1 ) / 2);
	vDSP_vclrD(target + center, 1, length - center);
}
C74_HIDDEN static t_class const * class = NULL;
C74_HIDDEN void*new(t_symbol const * const symbol, long const argc, t_atom const * const argv) {
	register t_pitcheffect * const this = (t_pitcheffect * const)object_alloc((t_class*const)class);
    if ( this ) {
		// init
		this->forward = vDSP_DFT_zop_CreateSetupD(NULL, 0, vDSP_DFT_FORWARD);
		this->inverse = vDSP_DFT_zop_CreateSetupD(NULL, 0, vDSP_DFT_INVERSE);
		this->cache = (double*const)sysmem_newptrclear(0);
		*(double*const)&this->pitch = 1;
		*(double*const)&this->width = 10;
		*(double*const)&this->verge = 1e-12;
		*(double*const)&this->level = 0;
		// attr
		attr_args_process(this, argc, (t_atom*const)argv);
		
		// init
        z_dsp_setup(&this->super, 2);
        outlet_new(this, "signal");
    }
    return this;
}
C74_HIDDEN void del(t_pitcheffect * const this) {
    z_dsp_free(&this->super);
	vDSP_DFT_DestroySetupD(this->inverse);
	vDSP_DFT_DestroySetupD(this->forward);
	sysmem_freeptr(this->cache);
}
C74_HIDDEN void dsp(t_pitcheffect * const this, double * const o, double const * const x, double const * const y, intptr_t const query, double const fs) {
	register intptr_t const count = this->count;
	register intptr_t const frame = this->frame;
	register intptr_t const index = this->index;
	register double const * const hanning = this->cache + frame * 0;
	register double const * const hilbert = this->cache + frame * 1;
	register double * const qr = this->cache + 2 * frame;
	register double * const qi = this->cache + 3 * frame;
	register double * const pr = this->cache + 4 * frame;
	register double * const pi = this->cache + 5 * frame;
	register double * const enque = this->cache + frame * 6;
	register double * const deque = this->cache + frame * 6 + count;
	register vDSP_DFT_SetupD const forward = this->forward;
	register vDSP_DFT_SetupD const inverse = this->inverse;
	{/* INPUT */
		register intptr_t const point = ( index + frame ) % count;
		register intptr_t const limit = ( count - point );
		register intptr_t const head = simd_reduce_min((simd_long2 const){0 + limit, query});
		register intptr_t const tail = simd_reduce_max((simd_long2 const){0, query - limit});
		sysmem_copyptr(y, enque + point, sizeof(double const) * head);
		sysmem_copyptr(y + limit, enque, sizeof(double const) * tail);
	}
	{/* PROCESS */
		register intptr_t phase = this->phase;
		for ( ; phase < query ; phase += fs / simd_reduce_max((simd_double2 const){1, fabs(x[phase])}) ) {
//			register intptr_t const width = simd_reduce_min((simd_long2 const){frame, 2 * fs / x[phase]});
			register intptr_t const point = ( index + phase ) % count;
			register intptr_t const limit = ( count ) - point;
			register intptr_t const head = simd_reduce_min((simd_long2 const){0 + limit, frame});
			register intptr_t const tail = simd_reduce_max((simd_long2 const){0, frame - limit});
			sysmem_copyptr(enque + point, qr, sizeof(double const) * head);
			sysmem_copyptr(enque, qr + limit, sizeof(double const) * tail);
			
			// initialize-windows
//			vDSP_vclrD(hanning, 1, frame);
//			vDSP_hann_windowD(hanning, width, vDSP_HANN_DENORM);
//			vDSP_vclrD(hilbert, 1, frame);
//			sysmem_copyptr(hanning + (width + 1) / 2, hilbert, sizeof(double const) * (width / 2 + 1));
//			vDSP_vsmulD(hilbert + 1, 1, (double const[]){2}, hilbert + 1, 1, (frame - 1 ) / 2);
			
			// hanning(t) * x(t) -> q(t)
			vDSP_vmulD(qr, 1, hanning, 1, qr, 1, frame);
			vDSP_vclrD(qi, 1, frame);
			
			// f(t) -> F(ω)
			vDSP_DFT_ExecuteD(forward, qr, qi, pr, pi);
			
			// F(ω) -> ln|F(ω)|
			vDSP_vmmaD(pr, 1, pr, 1, pi, 1, pi, 1, pr, 1, frame);
			vDSP_vthrD(pr, 1, &this->verge, pr, 1, frame);
			vvlog(pr, pr, (const int[]){(int const)frame});
			vDSP_vsdivD(pr, 1, (double const[]){2}, pr, 1, frame);
			vDSP_vclrD(pi, 1, frame);
			
			// ln|F(ω)| -> f(t)
			vDSP_DFT_ExecuteD(inverse, pr, pi, qr, qi);
			vDSP_vsdivD(qr, 1, (double const[]){(double const)frame}, qr, 1, frame);
//			vDSP_vsdivD(qi, 1, (double const[]){(double const)frame}, qi, 1, frame); // not necessary
			vDSP_vclrD(qi, 1, frame);
			
			// hilbert(t) * q(t) -> q(t)
			vDSP_vmulD(qr, 1, hilbert, 1, qr, 1, frame);
			
			// f(t) -> lnF(ω)
			vDSP_DFT_ExecuteD(forward, qr, qi, pr, pi);
			
			// lnF(ω) -> F(ω)
			vvexp(pr, pr, (const int[]){(int const)frame});
			vvsincos(pi, qi, qr, (const int[]){(int const)frame});
			vDSP_vmulD(qi, 1, pr, 1, pi, 1, frame);
			vDSP_vmulD(qr, 1, pr, 1, pr, 1, frame);
			
			// F(ω) -> f(t)
			vDSP_DFT_ExecuteD(inverse, pr, pi, qr, qi);
			vDSP_vsdivD(qr, 1, (double const[]){(double const)frame}, qr, 1, frame);
//			vDSP_vsdivD(qi, 1, (dou ble const[]){(double const)frame}, qi, 1, frame); // not necessary
//			vDSP_vclrD(pi, 1, frame); // not necessary
			
			vDSP_vaddD(qr, 1, deque + point, 1, deque + point, 1, head);
			vDSP_vaddD(qr + head, 1, deque, 1, deque, 1, tail);
		}
		this->phase = phase - query;
	}
	{/* OUTPUT */
		register intptr_t const point = index % count;
		register intptr_t const limit = count - point;
		register intptr_t const head = simd_reduce_min((simd_long2 const){0 + limit, query});
		register intptr_t const tail = simd_reduce_max((simd_long2 const){0, query - limit});
		sysmem_copyptr(deque + point, o, sizeof(double const) * head);vDSP_vclrD(deque + point, 1, head);
		sysmem_copyptr(deque, o + limit, sizeof(double const) * tail);vDSP_vclrD(deque, 1, tail);
	}
	{/* PROGRESS */
		this->index = ( index + query ) % count;
	}
}
C74_HIDDEN void nop(t_pitcheffect * const this, t_object const * const dsp64, double const * const * const ins, long const numins, double * const * const outs, long const numouts, long const framecount, long const flags, void const * const userparam) {
	vDSP_vclrD(outs[0], 1, framecount);
}
C74_HIDDEN void fix(t_pitcheffect * const this, t_object const * const dsp64, double const * const * const ins, long const numins, double * const * const outs, long const numouts, long const framecount, long const flags, void const * const userparam) {
	double * const buffer = (double*const)alloca(framecount * sizeof(double const));
	vDSP_vfillD(&this->pitch, buffer, 1, framecount);
	dsp(this, outs[0], buffer, ins[1], framecount, (double const)(uintptr_t const)userparam);
}
C74_HIDDEN void dyn(t_pitcheffect * const this, t_object const * const dsp64, double const * const * const ins, long const numins, double * const * const outs, long const numouts, long const framecount, long const flags, void const * const userparam) {
	dsp(this, outs[0], ins[0], ins[1], framecount, (double const)(uintptr_t const)userparam);
}
C74_HIDDEN void dsp64(t_pitcheffect * const this, t_object const * const dsp64, short const * const count, double const samplerate, long const maxvectorsize, long const flags) {
	register intptr_t const width = simd_reduce_max((simd_long2 const){1, samplerate * this->width / 1000});
	register intptr_t const frame = 2 << ilogb(width - 1);
	// initialize
	this->index = 0;
	this->frame = frame;
	this->phase = frame;
	this->count = simd_reduce_max((simd_ulong2 const){frame, maxvectorsize}) * 2;
	// resize fftsetup
	vDSP_DFT_DestroySetupD(this->forward);
	vDSP_DFT_DestroySetupD(this->inverse);
	this->forward = vDSP_DFT_zop_CreateSetupD(NULL, this->frame, vDSP_DFT_FORWARD);
	this->inverse = vDSP_DFT_zop_CreateSetupD(NULL, this->frame, vDSP_DFT_INVERSE);
	// memory allocation
	this->cache = (double*const)sysmem_resizeptrclear(this->cache, sizeof(double const) * (this->frame * 6 + this->count * 2));
	vDSP_vclrD(this->cache, 1, this->frame * 8 + this->count * 2);
	// hanning-window
	vDSP_hann_windowD(this->cache, width, vDSP_HANN_NORM);
	// hilbert-window
	hilbert(this->cache + this->frame, this->level, this->frame);
//	sysmem_copyptr(this->cache + (width + 1) / 2, this->cache + frame, sizeof(double const) * (width / 2 + 1));
//	vDSP_vsmulD(this->cache + frame + 1, 1, (double const[]){2}, this->cache + frame + 1, 1, (frame - 1 ) / 2);
	// entry
	dsp_add64((t_object*const)dsp64, (t_object*const)this, (t_perfroutine64 const)(count[1] ? count[0] ? dyn : fix : nop), 0, (uintptr_t const)samplerate);
}
C74_EXPORT void ext_main(void * const _) {
	C74_STATIC_ASSERT(sizeof(double const) <= sizeof(uintptr_t const), "pitchsense~ requires 64-bit architecture");
	if ( !class ) {
        t_class * const obj = class_new("pitcheffect~", (method const)new, (method const)del, sizeof(t_pitcheffect const), NULL, A_GIMME, 0);
		class_addattr(obj, attr_offset_new("pitch", gensym("float64"), 0, (method const)0L, (method const)0L, offsetof(t_pitcheffect, pitch)));
		class_addattr(obj, attr_offset_new("width", gensym("float64"), 0, (method const)0L, (method const)0L, offsetof(t_pitcheffect, width)));
		class_addattr(obj, attr_offset_new("threshold", gensym("float64"), 0, (method const)0L, (method const)0L, offsetof(t_pitcheffect, verge)));
		class_addattr(obj, attr_offset_new("pitchsuppress", gensym("float64"), 0, (method const)0L, (method const)0L, offsetof(t_pitcheffect, level)));
        class_addmethod(obj, (method const)dsp64, "dsp64", A_CANT, 0);
        class_dspinit(obj);
        class_register(CLASS_BOX, obj);
		class = obj;
    }
}
