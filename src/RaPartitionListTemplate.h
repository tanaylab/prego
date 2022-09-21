#ifndef stdalg_RaPartitionList_t
#define stdalg_RaPartitionList_t 1

#include "RaPartitionList.h"

template<class T>
void RaPartitionList<T>::remove_subset(uint subset_id)
{
	PartitionData &pd = m_subset_limits[subset_id];

	if(pd.size == 0) {
		return;
	}

	int cur = pd.get_front_id();
	int back_id = pd.get_back_id();
	while(cur != back_id) {
		m_list.remove(cur);
	}
	pd.size = 0;
}

template<class T>
void RaPartitionList<T>::remove_abs(int elem_id)
{
	int sbs_id = m_id_subset[elem_id];

	PartitionData &pd = m_subset_limits[sbs_id];
	if(pd.size == 1) {
		pd.size = 0;
	} else if(pd.get_front_id() == elem_id) {
		pd.set_front_id(m_list.get_next_id(elem_id));
	} else if(pd.get_back_id() == elem_id) {
		pd.set_back_id(m_list.get_prev_id(elem_id));
	}

	m_list.remove(elem_id);
}

template<class T>
int RaPartitionList<T>::pop_subset(uint subset_id)
{
	if(m_subset_limits.size() <= subset_id) {
		return(-1);
	}
	PartitionData &pd = m_subset_limits[subset_id];
	if(pd.size == 0) {
		return(-1);
	}
	int last = pd.get_back_id();
	if(pd.size == 1) {
		pd.set_back_id(-1);
		pd.set_front_id(-1);
	} else {
		pd.set_back_id(m_list.get_prev_id(last));
	}
	pd.size--;
	return(last);
}

//allocating id if needed
template<class T>
int RaPartitionList<T>::push_back_subset(uint subset_id, const T &val, int elem_id)
{
	if(m_subset_limits.size() <= subset_id) {
		m_subset_limits.resize(subset_id + 1);
	}
	PartitionData &pd = m_subset_limits[subset_id];
	if(pd.size == 0) {
		elem_id = m_list.push_back(val, elem_id);
		if(m_id_subset.size() <= uint(elem_id)) {
			m_id_subset.resize(elem_id + 1);
		}
		m_id_subset[elem_id] = subset_id;
		pd.set_front_id(elem_id);
	} else {
		elem_id = m_list.insert_after(subset_id, pd.get_back_id(), elem_id);
		if(m_id_subset.size() <= uint(elem_id)) {
			m_id_subset.resize(elem_id + 1);
		}
		m_id_subset[elem_id] = subset_id;
	}
	pd.set_back_id(elem_id);
	pd.size++;
	return(elem_id);
}

#endif // stdalg_ds_RaPartitionList_t
