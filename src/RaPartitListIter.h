#ifndef stdalg_ds_RaPartitListIter_h
#define stdalg_ds_RaPartitListIter_h 1

#include "RaPartitionList.h"

template<class T>
class RaPartitListIter {

protected:

	const RaPartitionList<T> &m_list;

	int m_cur_subset;

	uint m_last_id;

	int m_cur_id;

public:

	RaPartitListIter(const RaPartitionList<T> &list, int subset_id = -1) :
		m_list(list),
		m_cur_subset(subset_id)
	{
		if(subset_id != -1) {
			first();
		}
	}
	
	void first() {
		const typename RaPartitionList<T>::PartitionData &pd = 
				m_list.m_subset_limits[m_cur_subset];
		if(pd.size == 0) {
			m_cur_id = -1;
		} else {
			m_cur_id  = pd.get_front_id();
			m_last_id = pd.get_back_id();
		}
	}

	void restart(int subset_id) {
		m_cur_subset = subset_id;
		first();
	}

	void operator++() {
		if(m_cur_id == int(m_last_id)) {
			m_cur_id = -1;
		} else {
			m_cur_id = m_list.get_next_id(m_cur_id);
		}
	}

	bool is_valid() const {
		return(m_cur_id != -1);
	}

	int get_cur_id() {
		return(m_cur_id);
	}
};

#endif // stdalg_ds_RaPartitListIter_h
