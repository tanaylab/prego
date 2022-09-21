#ifndef stdalg_ds_RaList_t
#define stdalg_ds_RaList_t 1

#include "RaList.h"
#include "RaListIter.h"

template<class T> 
int RaList<T>::push_front(const T &obj, int elem_id)
{
	if(elem_id == -1) {
		elem_id = get_free_id();
	}
	if(uint(elem_id) >= m_nodes.size()) {
		m_nodes.resize(elem_id + 1);
	}
	RaListNode<T> &node = m_nodes[elem_id];
	node.get_obj() = obj;
	node.id = elem_id;
	node.prev_id = -1;
	node.next_id = m_front_id;
	if(m_size == 0) {
		m_back_id = elem_id;
	} else {
		m_nodes[m_front_id].prev_id = elem_id;
	}
	m_front_id = elem_id;
	m_size++;
	return(elem_id);
}

template<class T> 
int RaList<T>::push_back(const T &obj, int elem_id)
{
	if(elem_id == -1) {
		elem_id = get_free_id();
	}
	if(uint(elem_id) >= m_nodes.size()) {
		m_nodes.resize(elem_id + 1);
	}
	RaListNode<T> &node = m_nodes[elem_id];
	node.get_obj() = obj;
	node.id = elem_id;
	node.prev_id = m_back_id;
	node.next_id = -1;
	if(m_size == 0) {
		m_front_id = elem_id;
	} else {
		m_nodes[m_back_id].next_id = elem_id;
	}
	m_back_id = elem_id;
	m_size++;
	return(elem_id);
}

template<class T> 
int RaList<T>::insert_after(const T &obj, int after_id, int elem_id)
{
	if(elem_id == -1) {
		elem_id = get_free_id();
	}
	if(uint(elem_id) >= m_nodes.size()) {
		m_nodes.resize(elem_id + 1);
	}
	RaListNode<T> &node = m_nodes[elem_id];
	RaListNode<T> &after = m_nodes[after_id];
	node.get_obj() = obj;
	node.id = elem_id;
	node.prev_id = after_id;
	node.next_id = after.next_id;
	after.next_id = elem_id;

	if(after_id == m_back_id) {
		m_back_id = elem_id;
	}
	m_size++;
	return(elem_id);
}

template<class T> 
void RaList<T>::remove(int elem_id)
{
	m_size--;
	if(m_size == 0) {
		m_front_id = -1;
		m_back_id = -1;
		return;
	}
	RaListNode<T> &node = m_nodes[elem_id];
	node.id = -1;

	if(m_front_id == elem_id) {
		m_front_id = node.next_id;
		m_nodes[m_front_id].prev_id = -1;
	} else if(m_back_id == elem_id) {
		m_back_id = node.prev_id;
		m_nodes[m_back_id].next_id = -1;
	} else {
		m_nodes[node.prev_id].next_id = node.next_id;
		m_nodes[node.next_id].prev_id = node.prev_id;
	}
}

template<class T>
int RaList<T>::get_free_id()
{
	//the easy way out.., we should handle a list of free ids
	return(m_nodes.size());
}

template<class T>
RaListIter<T> RaList<T>::begin()
{
	return(RaListIter<T>(*this));
}

template<class T>
RaListIter<T> RaList<T>::end() 
{
	RaListIter<T> i(*this);
	i.last();
	return(i);
}

#endif // stdalg_ds_RaList_t
