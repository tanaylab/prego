#ifndef stdalg_ds_RaPartitionList_h
#define stdalg_ds_RaPartitionList_h 1

#include "RaList.h"

template<class T> class RaPartitListIter;

template<class T>
class RaPartitionList {

public:

	friend class RaPartitListIter<T>;

private:

	class PartitionData {
	
	private:
		int m_front_id;
		int m_back_id;
	public:
		int size;

		int get_front_id() const {
			return(m_front_id);
		}
		int set_front_id(int f) {
			return(m_front_id = f);
		}
		int get_back_id() const {
			return(m_back_id);
		}
		int set_back_id(int b) {
			return(m_back_id = b);
		}

		PartitionData() :
			size(0)
		{}
	};

public:
	typedef RaPartitListIter<T> iterator;

protected:

	RaList<T> m_list;

	std::vector<PartitionData> m_subset_limits;
	std::vector<uint> m_id_subset;

public:

	void remove_subset(uint subset_id);
	void remove_abs(int elem_id);

	bool empty_susbset(uint subset_id) const {
		return(subset_id >= m_subset_limits.size()
		    || m_subset_limits[subset_id].size == 0);
	}

	int get_subset_head(uint subset_id) const {
		return(subset_id < m_subset_limits.size() ? 
		    m_subset_limits[subset_id].get_front_id() : -1);
	}

	int pop_subset(uint subset_id);

	//allocating id if needed
	int push_back_subset(uint subset_id, const T &val, int elem_id = -1);

	const T &operator[](int elem_id) {
		return(m_list[elem_id]);
	}

	int get_next_id(int elem_id) const {
		return(m_list.get_next_id(elem_id));
	}
	int get_prev_id(int elem_id) const {
		return(m_list.get_prev_id(elem_id));
	}

	bool is_member(int elem_id) const {
		return(m_list.is_member(elem_id));
	}
};

#endif // stdalg_ds_RaPartitionList_h
