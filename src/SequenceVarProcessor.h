#ifndef SEQUENCEVARPROCESSOR_H_
#define SEQUENCEVARPROCESSOR_H_

#include "TrackExpressionVars.h"
#include "GInterval.h"
#include "GenomeSeqFetch.h"
#include "rdbutils.h"
#include <vector>

class SequenceVarProcessor {
public:
	SequenceVarProcessor(rdb::IntervUtils &iu, GenomeSeqFetch &shared_seqfetch)
		: m_iu(iu), m_shared_seqfetch(shared_seqfetch) {}

	void process_sequence_vars(
		TrackExpressionVars::Track_vars &track_vars,
		const GInterval &interval,
		unsigned idx
	);

private:
	rdb::IntervUtils &m_iu;
	GenomeSeqFetch &m_shared_seqfetch;

	void batch_process_sequence_vtracks(
		std::vector<TrackExpressionVars::Track_var*> &kmer_vtracks,
		std::vector<TrackExpressionVars::Track_var*> &pwm_vtracks,
		const GInterval &interval,
		unsigned idx
	);

	void process_individual_sequence_vars(
		std::vector<TrackExpressionVars::Track_var*> &kmer_vtracks,
		std::vector<TrackExpressionVars::Track_var*> &pwm_vtracks,
		std::vector<TrackExpressionVars::Track_var*> &masked_vtracks,
		const GInterval &interval,
		unsigned idx
	);
};

#endif /* SEQUENCEVARPROCESSOR_H_ */

