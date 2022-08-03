#ifndef stdalg_ds_RaListIter_h
#define stdalg_ds_RaListIter_h 1

#include <vector>
#include "RaList.h"

template<class T> 
class RaListIter {

protected:

	const RaList<T> &m_ralist;

	const RaListNode<T> *cur_node;

	int m_cur_id;

public:

	RaListIter(RaList<T> &papa) :
		m_ralist(papa)
	{}

	bool is_valid() const {
		return(m_cur_id != -1);
	}

	void operator++() {
		m_cur_id = cur_node->next_id;
		cur_node = &(m_ralist.m_nodes[m_cur_id]);
	}
	void operator--() {
		m_cur_id = cur_node->prev_id;
		cur_node = &(m_ralist.m_nodes[m_cur_id]);
	}

	const RaListNode<T> &operator->() {
		return(*cur_node);
	}

	const RaListNode<T> &operator*() {
		return(*cur_node);
	}

	int get_id() const {
		return(m_cur_id);
	}

	void first() {
		m_cur_id = m_ralist.get_front_id();
	}
	void last() {
		m_cur_id = m_ralist.get_back_id();
	}
};

template<class T>
bool operator==(RaListIter<T> &o1, RaListIter<T> &o2)
{
	return(o1.get_id() == o2.get_id());
}

#endif // stdalg_ds_RaPartitionList_h
